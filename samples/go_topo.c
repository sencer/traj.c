#include "xyz.h"
#include "topo.h"
#include "periodic_table.h"

void ConvertAtomIDsToTypes(Crystal *c)
{
  for (int i = 0; i < c->nat; ++i)
  {
    switch (c->atoms[i].Z)
    {
      case 1: // Hydrogen
        switch (c->atoms[i].id)
        {
          case 1:
            c->atoms[i].id = 8; break;
          case 2:
            c->atoms[i].id = 9; break;
          case 3:
            c->atoms[i].id = 10; break;
          default:
            fprintf(stderr, "H%d is undefined.\n", c->atoms[i].id);
        }
        break;

      case 6: // Carbon
        switch (c->atoms[i].id)
        {
          case 1:
            c->atoms[i].id = 0; break;
          case 2:
            c->atoms[i].id = 1; break;
          case 3:
            c->atoms[i].id = 2; break;
          case 4:
            c->atoms[i].id = 3; break;
          default:
            fprintf(stderr, "C%d is undefined.\n", c->atoms[i].id);
        }
        break;

      case 8: // Oxygen
        switch (c->atoms[i].id)
        {
          case 1:
            c->atoms[i].id = 4; break;
          case 2:
            c->atoms[i].id = 5; break;
          case 3:
            c->atoms[i].id = 6; break;
          case 4:
            c->atoms[i].id = 7; break;
          default:
            fprintf(stderr, "O%d is undefined.\n", c->atoms[i].id);
        }
        break;

      default:
        fprintf(stderr, "%s is undefined.\n", PT_Symbol(c->atoms[i].Z));
    }
  }
}

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
  return (dist<1.2||(typ>2 && dist<1.4)||(typ>10 && dist<1.95))?1:0;
}

BondingInfo *BondSearch(Crystal *c)
{
  CoarseBox *box = BoxInit(c, 2.0);
  BoxFill(c, box);
  BondingInfo *bnd = BondingInit(c);
  BondingPopulate(c, box, bnd, CheckBonding);
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
  // Will read a GO structure as an xyz file and generate some topology
  // information.
  // The second line of the xyz should include "celldm x y z"
  // Atoms in the xyz file should be named C1-C4,O1-O4,H1-H2 as defined below.
  // First column is the index corresponding to each in this source code.
  // 0 C1 -> sp2 carbon          CA
  // 1 C2 -> C-oh                CF
  // 2 C3 -> epoxy carbon        CT
  // 3 C4 -> carboxyl carbon     C
  // 4 O1 -> c-Oh                OH
  // 5 O2 -> epoxy O             OS
  // 6 O3 -> carboxyl O          O_3
  // 7 O4 -> carboxyl Oh         OH2
  // 8 H1 -> c-oH                HO
  // 9 H2 -> carboxyl H          HO2
  //10 H3 -> benzene  H          H

  // we will use up to four atoms (when writing Dihedrals for example)
  // so lets allocate some memory for them.
  int typ1, typ2, typ3, typ4,
      // some memory for atm ids also
      atm1, atm2, atm3, atm4;

  char
    // we need to teach opls atom types for each of these atoms
    opls_type[11][9] = {
      "opls_147", "opls_166", "opls_182", "opls_267", "opls_167", "opls_179",
      "opls_269", "opls_268", "opls_168", "opls_270", "opls_146"
    },
    // we will also want to output the atom types in .top file for human reader
    defs[11][4] = {
      "CA", "CF", "CT", "C", "OH", "OS", "O", "OH2", "HO", "HO2", "H"
    };

  double
    // charges of each atom type, taken from literature TODO where exactly?
    charges[11] = {
      0.0,  0.15,  0.14,  0.52, -0.585, -0.28, -0.44, -0.53, 0.435, 0.45, 0.00
    },
    // let's keep the total charge to make sure the whole molecule is neutral
    tot_charge,
    // masses of each atom type
    masses[11] = {
      12.01100, 12.01100, 12.01100, 12.01100, 15.99940, 15.99940, 15.99940,
      15.99940, 1.00794,  1.00794,  1.00794
    };

  // Read the xyz file to a Crystal
  Crystal *c = ReadXYZ(argv[1]);
  ConvertAtomIDsToTypes(c);

  // Open the argv[2] to write the .top file
  FILE *top = fopen(argv[2], "w");

  fprintf(top, "#include \"oplsaa.ff/forcefield.itp\"\n");
  fprintf(top, "\n[ moleculetype ]\n");
  fprintf(top, "; name  nrexcl\n");
  fprintf(top, "     GO   3\n");

  // Now let's write the [atoms] section
  fprintf(top, "\n[ atoms ]\n");
  for (int i = 0; i < c->nat; ++i)
  {
    typ1 = c->atoms[i].id;
    tot_charge += charges[typ1];
    fprintf(top, "%5d    %-10s   1     GO_s      %s %4d %9.5f %9.5f ; %9.5f\n",
        i+1, opls_type[typ1], PT_Symbol(c->atoms[i].Z), i+1, charges[typ1],
        masses[typ1], tot_charge);
  }

  // for all other sections, we will need to know the bonding structure
  BondingInfo *bnd = BondSearch(c);

  // Print out the bonds section
  fprintf(top, "\n[ bonds ]\n");
  // iterate through each atom -> atm1
  for (atm1 = 0; atm1 < c->nat; ++atm1)
  {
    typ1 = c->atoms[atm1].id;
    // iterate through each bond atm1 makes -> atm2
    for (int i = 0; i < bnd->nbonds[atm1]; ++i)
    {
      atm2 = bnd->bonds[atm1][i];
      // Look at each bond only once
      if (atm1 < atm2)
      {
        typ2 = c->atoms[atm2].id;
        fprintf(top, "   %5d %5d   1                        ; %s-%s\n", 1+atm1,
            1+atm2, defs[typ1], defs[typ2]);
      }
    }
  }

  // calculate the pairs using bonding information
  // data will be hold in a matrix where (i,j)th element gives the minimum
  // number of bonds needed to traverse to go from i to j.
  // since this is a symmetric matrix, we will keep only half of it
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
  for (int k = 0; k < c->nat; ++k)
  {
    for (int i = 0; i < c->nat; ++i)
    {
      ik = I(c->nat, i, k);
      for (int j = i+1; j < c->nat; ++j)
      {
        ij = I(c->nat, i, j);
        kj = I(c->nat, k, j);
        if (pairs[ij] > pairs[ik]+pairs[kj])
        {
          pairs[ij] = pairs[ik] + pairs[kj];
        }
      }
    }
  }

  // and print the 1-4 pairs only
  fprintf(top, "\n[ pairs ]\n");
  for (int atm1 = 0; atm1 < c->nat; ++atm1)
  {
    for (int atm2 = atm1+1; atm2 < c->nat; ++atm2)
    {
      if (pairs[I(c->nat, atm1, atm2)] == 3)
      {
        fprintf(top, "   %5d %5d   1; %s-%s\n", 1+atm1, 1+atm2,
            defs[c->atoms[atm1].id], defs[c->atoms[atm2].id]);
      }
    }
  }
  free(pairs);

  // Print out angles section
  fprintf(top, "\n[ angles ]\n");

  int nangles;
  int *angles = Angles(c, bnd, &nangles);

  for (int i = 0; i < nangles; ++i)
  {
    typ1 = c->atoms[angles[3*i]].id;
    typ2 = c->atoms[angles[3*i+1]].id;
    typ3 = c->atoms[angles[3*i+2]].id;
    fprintf(top, "%4d %4d %4d 1                       ; %s-%s-%s\n",
        1+angles[3*i], 1+angles[3*i+1], 1+angles[3*i+2], defs[typ1],
        defs[typ2], defs[typ3]);
  }
  free(angles);

  // Print out dihedrals section
  fprintf(top, "\n[ dihedrals ]\n");

  int ndiheds;
  int *diheds = Dihedrals(c, bnd, &ndiheds);

  for (int i = 0; i < ndiheds; ++i)
  {
    typ1 = c->atoms[diheds[4*i]].id;
    typ2 = c->atoms[diheds[4*i+1]].id;
    typ3 = c->atoms[diheds[4*i+2]].id;
    typ4 = c->atoms[diheds[4*i+3]].id;
    fprintf(top, "%5d %5d %5d %5d 3           ; %s-%s-%s-%s\n",
        1+diheds[4*i], 1+diheds[4*i+1], 1+diheds[4*i+2], 1+diheds[4*i+3],
        defs[typ4], defs[typ1], defs[typ2], defs[typ3]);
  }
  free(diheds);

  // Print out impropers section
  fprintf(top, "\n[ dihedrals ] ; impropers\n");

  int nimprp;
  int *imprps = Impropers(c, bnd, &ndiheds);

  for (int i = 0; i < ndiheds; ++i)
  {
    typ1 = c->atoms[imprps[4*i]].id;
    typ2 = c->atoms[imprps[4*i+1]].id;
    typ3 = c->atoms[imprps[4*i+2]].id;
    typ4 = c->atoms[imprps[4*i+3]].id;
    fprintf(top, "%5d %5d %5d %5d 3             ; %s-%s-%s-%s\n",
        1+imprps[4*i], 1+imprps[4*i+1], 1+imprps[4*i+2], 1+imprps[4*i+3],
        defs[typ4], defs[typ1], defs[typ2], defs[typ3]);
  }
  free(imprps);

  fprintf(top, "\n[ system ]\n");
  fprintf(top, "; title from mol2 input\n");
  fprintf(top, "GO\n");
  fprintf(top, "\n[ molecules ]\n");
  fprintf(top, "; molecule name    nr.\n");
  fprintf(top, "     GO           1\n");

  BondingDelete(bnd);
  CrystalDelete(c);
  fclose(top);
  return 0;
}

// vi: fdm=syntax
