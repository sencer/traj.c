#include "xyz.h"
#include "topo.h"
#include "periodic_table.h"
#include "unistd.h"

Crystal *ReadXYZ(char *filename)
{
  FILE *xyz = fopen(filename, "r");
  int nat;
  double dm[3];
  XYZReadHeader(xyz, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  XYZReadFrame(xyz, c);
  fclose(xyz);
  return c;
}

int CheckBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.2||(typ>2 && dist<1.4)||(typ>10 && dist<1.84))?1:0;
}

BondingInfo *BondSearch(Crystal *c)
{
  CoarseBox *box = BoxInit(c, 2.0);
  BoxFill(c, box);
  BondingInfo *bnd = BondingInit(c);
  BondingPopulate(c, box, bnd, CheckBonding, 1);
  BoxDelete(box);
  return bnd;
}

int I(int nat, int i, int j)
{
  // a function to calculate the index of (i, j)th element of half diagonal
  // matrix in the linear memory. (i, j) == (j, i)
  return ((i<j)?(i*nat - i*(i+1)/2 + j):(j*nat - j*(j+1)/2 + i));
}


int main(int argc, char *argv[])
{

  // parse arguments
  char *fxyz, *ftop, *fbonds;
  int opt, ibonds = 0, idenoted = 0;
  while ((opt = getopt(argc, argv, "di:o:b:")) != -1)
  {
    switch (opt)
    {
      case 'i':
        fxyz = optarg;
        break;
      case 'o':
        ftop = optarg;
        break;
      case 'b':
        fbonds = optarg;
        ibonds = 1;
        break;
      case 'd':
        idenoted = 1;
        break;
      default:
        fprintf(stderr, "Usage: %s -i xyz_file -o top_file [-b bonds]\n", argv[0]);
        exit(-1);
    }
  }

  int atm1, atm2, atm3, atm4;

  // Read the xyz file to a Crystal
  Crystal *c = ReadXYZ(fxyz);
  // Determine or Read the chemical bonding bonding structure
  BondingInfo *bnd;
  if (ibonds)
  {
    bnd = BondingInit(c);
    BondingReadFile(fbonds, bnd);
    fprintf(stderr, "\rBonds read from file.%30s", "");
  }
  else
  {
    bnd = BondSearch(c);
    fprintf(stderr, "\rBonds calculated from distance search.%30s", "");
  }

  // Open the argv[2] to write the .top file
  FILE *top = fopen(ftop, "w");

  fprintf(top, "#include \"oplsaa.ff/forcefield.itp\"\n");
  fprintf(top, "\n[ moleculetype ]\n");
  fprintf(top, "; name  nrexcl\n");
  fprintf(top, "     Tz4   3\n");

  // Print out the bonds section
  fprintf(top, "\n[ bonds ]\n");
  // iterate through each atom -> atm1
  for (atm1 = 0; atm1 < c->nat; ++atm1)
  {
    // iterate through each bond atm1 makes -> atm2
    for (int i = 0; i < bnd->nbonds[atm1]; ++i)
    {
      atm2 = bnd->bonds[atm1][i];
      // Look at each bond only once
      if (atm1 < atm2)
      {
        fprintf(top, "   %5d %5d   1                        ; %s-%s\n", 1+atm1,
            1+atm2, PT_Symbol(c->atoms[atm1].Z), PT_Symbol(c->atoms[atm2].Z));
      }
    }
  }

  fprintf(stderr, "\rBonds written into topology file.%30s", "");
  // calculate the pairs using bonding information
  // data will be hold in a matrix where (i,j)th element gives the minimum
  // number of bonds needed to traverse to go from i to j.
  // since this is a symmetric matrix, we will keep only half of it
  fprintf(stderr, "\rCalculating pairs.%30s", "");
  int *pairs = malloc(c->nat*(1+c->nat)/2*sizeof(int));
  int ij, ik, kj;

  // Matrix will be filled using Floyd-Warshall
  // Initialize all distances to "infinity" except for diagonals, which are 0
  for (int row = 0; row < c->nat; ++row)
  {
    // start from the diagonal, this is its index:
    ij = I(c->nat, row, row);
    pairs[ij] = 0;
    // then fill rest of the row
    for (int i = 1; i < c->nat-row; ++i)
    {
      pairs[ij+i] = 10000;
    }
  }

  // Fill in the bonding information
  for (atm1 = 0; atm1 < c->nat; ++atm1)
  {
    for (int i = 0; i < bnd->nbonds[atm1]; ++i)
    {
      atm2 = bnd->bonds[atm1][i];
      if (atm1<atm2)
      {
        pairs[I(c->nat, atm1, atm2)] = 1;
      }
    }
  }

  // Now calculate the real distances
  // this is the heaviest part of the calculation
  // in practice, only 1-4 bonds are needed, and CoarseBox idea can be used
  // there too. Perhaps one day. TODO
#pragma omp parallel for private(ik,ij,kj)
  for (int k = 0; k < c->nat; ++k)
  {
    for (int i = 0; i < c->nat; ++i)
    {
      ik = I(c->nat, i, k);
      if (pairs[ik] == 10000||i==k) continue;
      for (int j = i+1; j < c->nat; ++j)
      {
        if (j==k) continue;
        ij = I(c->nat, i, j);
        kj = I(c->nat, k, j);
        if (pairs[kj]!=10000 && pairs[ij] > pairs[ik]+pairs[kj])
        {
          pairs[ij] = pairs[ik] + pairs[kj];
        }
      }
    }
  }

  fprintf(stderr, "\rPairs calculated.%30s", "");
  // and print the 1-4 pairs only
  fprintf(top, "\n[ pairs ]\n");
  for (int atm1 = 0; atm1 < c->nat; ++atm1)
  {
    for (int atm2 = atm1+1; atm2 < c->nat; ++atm2)
    {
      if (pairs[I(c->nat, atm1, atm2)] == 3)
      {
        fprintf(top, "   %5d %5d   1; %s-%s\n", 1+atm1, 1+atm2,
            PT_Symbol(c->atoms[atm1].Z), PT_Symbol(c->atoms[atm2].Z));
      }
    }
  }
  free(pairs);
  fprintf(stderr, "\rPairs written to topology file.%30s", "");

  // Print out angles section
  fprintf(top, "\n[ angles ]\n");

  int nangles;
  int *angles = Angles(c, bnd, &nangles);

  for (int i = 0; i < nangles; ++i)
  {
    atm1 = angles[3*i];
    atm2 = angles[3*i+1];
    atm3 = angles[3*i+2];

    fprintf(top, "%4d %4d %4d 1", 1+atm1, 1+atm2, 1+atm3);
    fprintf(top, "%23s", "");
    fprintf(top, "; %s-%s-%s\n", PT_Symbol(c->atoms[atm1].Z), PT_Symbol(c->atoms[atm2].Z), PT_Symbol(c->atoms[atm3].Z));
  }
  free(angles);
  fprintf(stderr, "\rAngles written to topology file.%30s", "");

  // Print out dihedrals section
  fprintf(top, "\n[ dihedrals ]\n");

  int ndiheds;
  int *diheds = Dihedrals(c, bnd, &ndiheds);

  for (int i = 0; i < ndiheds; ++i)
  {
    atm1 = diheds[4*i];
    atm2 = diheds[4*i+1];
    atm3 = diheds[4*i+2];
    atm4 = diheds[4*i+3];

    fprintf(top, "%5d %5d %5d %5d 1", 1+atm1, 1+atm2, 1+atm3, 1+atm4);
    fprintf(top, "%67s", "");
    fprintf(top, "; %s-%s-%s-%s\n", PT_Symbol(c->atoms[atm1].Z), PT_Symbol(c->atoms[atm2].Z), PT_Symbol(c->atoms[atm3].Z), PT_Symbol(c->atoms[atm4].Z));
  }
  free(diheds);
  fprintf(stderr, "\rDihedrals written to topology file.%30s", "");

  // Print out impropers section
  fprintf(top, "\n[ dihedrals ] ; impropers\n");

  int nimprp;
  int *imprps = Impropers(c, bnd, &ndiheds);

  for (int i = 0; i < ndiheds; ++i)
  {
    atm1 = imprps[4*i];
    atm2 = imprps[4*i+1];
    atm3 = imprps[4*i+2];
    atm4 = imprps[4*i+3];
    fprintf(top, "%5d %5d %5d %5d 2             ; %s-%s-%s-%s\n",
        1+atm1, 1+atm2, 1+atm3, 1+atm4,
        PT_Symbol(c->atoms[atm1].Z), PT_Symbol(c->atoms[atm2].Z), PT_Symbol(c->atoms[atm3].Z), PT_Symbol(c->atoms[atm4].Z));
  }
  free(imprps);
  fprintf(stderr, "\rImpropers written to topology file.%30s", "");

  fprintf(top, "\n[ system ]\n");
  fprintf(top, "; title from mol2 input\n");
  fprintf(top, "Tz4\n");
  fprintf(top, "\n[ molecules ]\n");
  fprintf(top, "; molecule name    nr.\n");
  fprintf(top, "     GO           1\n");

  BondingDelete(bnd);
  CrystalDelete(c);
  fclose(top);
  fprintf(stderr, "\rDONE%50s\n", "");
  return 0;
}

// vi: fdm=syntax
