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

int DihedralParams(int i, int j, int k, int l, double *params)
{
  int index = 729*i + 81*j + 9*k + l;
  if ( index == 162 || index == 18 || index == 23 || index == 3807)
  {
    params[0] = 0.00000;
    params[1] = 0.00000;
    params[2] = 0.00000;
    params[3] = 0.00000;
    params[4] = 0.00000;
    params[5] = 0.00000;
  }
  else if ( index == 1881 || index == 209)
  {
    params[0] = 1.71544;
    params[1] = 2.84512;
    params[2] = 1.04600;
    params[3] = -5.60656;
    params[4] = 0.00000;
    params[5] = 0.00000;
  }
  else if ( index == 125 || index == 6165)
  {
    params[0] =  7.03749;
    params[1] =  0.00000;
    params[2] =  -7.03749;
    params[3] =  0.00000;
    params[4] =  0.00000;
    params[5] =  0.00000;
  }
  else if ( index == 3159 || index == 31 || index == 33 || index == 34 || index
      == 4617 || index == 4618 || index == 5346 || index == 762)
  {
    params[0] = .78640;
    params[1] = .00000;
    params[2] = 8.78640;
    params[3] = .00000;
    params[4] = .00000;
    params[5] = .00000;
  }
  else if ( index == 4661 || index == 6189)
  {
    params[0] = 23.01200;
    params[1] = 0.00000;
    params[2] = -23.01200;
    params[3] = 0.00000;
    params[4] = 0.00000;
    params[5] = 0.00000;
  }
  else if ( index == 287 || index == 315 || index == 6183 || index == 7155)
  {
    params[0] =  29.28800;
    params[1] =  -8.36800;
    params[2] = -20.92000;
    params[3] =   0.00000;
    params[4] =   0.00000;
    params[5] =   0.00000;
  }
  else if ( index == 0 || index == 1 || index == 1459 || index == 1460 || index
      == 1467 || index == 2187 || index == 3 || index == 729 || index == 730 ||
      index == 731 || index == 738 || index == 81 || index == 82 || index == 83
      || index == 9)
  {
    params[0] = 30.33400;
    params[1] = 0.00000;
    params[2] = -30.33400;
    params[3] = 0.00000;
    params[4] = 0.00000;
    params[5] = 0.00000;
  }
  else if ( index == 1458 || index == 2)
  {
    params[0] =  44.97800;
    params[1] =  0.00000;
    params[2] =  -44.97800;
    params[3] =  0.00000;
    params[4] =  0.00000;
    params[5] =  0.00000;
  }
  else
  {
    params[0] =  -1;
    params[1] =  -1;
    params[2] =  -1;
    params[3] =  -1;
    params[4] =  -1;
    params[5] =  -1;
  }
  return 0;
}

double angle_tho(int i, int j, int k)
{
  int index = 81 * i + 9 * j + k;
  //debug
  if (
      index == 164 /*CT-CA-CT*/ || 
      index == 163 /*CT-CA-CF*/ ||
      index == 83  /*CF-CA-CT*/
     )
  {
    return 119.7;
  }
  else if (index == 82 /*CF-CA-CF*/)
  {
    return 120.0;
  }
  else
  {
    return -1;
  }
  // end debug: only write those not defined in GROMACS
  if ( index == 18  /* CA-CT-CA */ ||
       index == 23  /* CA-CT-OS */ ||
       index == 209 /* CT-OS-CT */ ||
       index == 423 /* OS-CT-CA */ )
  {
    return 109.5;
  }
  else if( index == 125 /*  */ ||
           index == 315 /*  */ ||
           index == 685 /*  */ ||
           index == 795 /*  */ )
  {
    return 113;
  }
  else if(index == 164)
  {
    return 116;
  }
  else if(index == 163 || index == 165 || index == 245 || index == 83)
  {
    return 119.7;
  }
  else if(index == 0 || index == 1 || index == 2 || index == 3 || index == 9 ||
      index == 13 || index == 31 || index == 34 || index == 81 || index == 82
      || index == 84 || index == 162 || index == 243 || index == 244 || index
      == 333 || index == 351 || index == 594)
  {
    return 120;
  }
  else if(index == 33 || index == 513)
  {
    return 120.4;
  }
  else if(index == 520 || index == 600)
  {
    return 121;
  }
  else
  {
    return -1;
  }
}

double angle_cth(int i, int j, int k)
{
  int index = 81 * i + 9 * j + k;
  //debug
  if(index==164||index == 163||index == 83)
  {
    return 585.76;
  }
  else if (index == 82)
  {
    return 527.184;
  }
  else
  {
    return -1;
  }

  if (index == 125 || index == 315 || index == 685 || index == 795 )
  {
    return 292.88;
  }
  else if (index == 18 )
  {
    return 334.72;
  }
  else if (index == 23 || index == 423 )
  {
    return 418.4;
  }
  else if (index == 209 )
  {
    return 502.08;
  }
  else if (index == 0 || index == 1 || index == 81 || index == 82 || index == 9 )
  {
    return 527.184;
  }
  else if (index == 13 || index == 162 || index == 163 || index == 164 || index
      == 165 || index == 245 || index == 2 || index == 31 || index == 333 ||
      index == 34 || index == 351 || index == 594 || index == 83 )
  {
    return 585.76;
  }
  else if (index == 33 || index == 513 || index == 520 || index == 600 )
  {
    return 669.44;
  }
  else if (index == 243 || index == 244 || index == 3 || index == 84 )
  {
    return 711.28;
  }
  else
  {
    return -1;
  }
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

  // a dummy variable to hold some force field parameters
  double params[6];

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
    atm1 = angles[3*i];
    atm2 = angles[3*i+1];
    atm3 = angles[3*i+2];
    typ1 = c->atoms[atm1].id;
    typ2 = c->atoms[atm2].id;
    typ3 = c->atoms[atm3].id;

    params[0] = angle_tho(typ1, typ2, typ3);

    if (params[0] == -1)
    {
      fprintf(top, "%4d %4d %4d 1", 1+atm1, 1+atm2, 1+atm3);
      fprintf(top, "%23s", "");
      fprintf(top, "; %s-%s-%s\n", defs[typ1], defs[typ2], defs[typ3]);
    }
    else
    {
      params[1] = angle_cth(typ1, typ2, typ3);
      fprintf(top, "%4d %4d %4d 1", 1+atm1, 1+atm2, 1+atm3);
      fprintf(top, " %10.5f %10.1f ", params[0], params[1]);
      fprintf(top, "; %s-%s-%s\n", defs[typ1], defs[typ2], defs[typ3]);
    }
  }
  free(angles);

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
    typ1 = c->atoms[atm1].id;
    typ2 = c->atoms[atm2].id;
    typ3 = c->atoms[atm3].id;
    typ4 = c->atoms[atm4].id;

    DihedralParams(typ1, typ2, typ3, typ4, params);
    if (params[0] == -1)
    {
      fprintf(top, "%5d %5d %5d %5d 1", 1+atm1, 1+atm2, 1+atm3, 1+atm4);
      fprintf(top, "%67s", "");
      fprintf(top, "; %s-%s-%s-%s\n", defs[typ1], defs[typ2], defs[typ3],
          defs[typ4]);
    }
    else
    {
      fprintf(top, "%5d %5d %5d %5d 3", 1+atm1, 1+atm2, 1+atm3, 1+atm4);
      fprintf(top, " %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f ", params[0],
          params[1], params[2], params[3], params[4], params[5]);
      fprintf(top, "; %s-%s-%s-%s\n", defs[typ1], defs[typ2], defs[typ3],
          defs[typ4]);
    }
  }
  free(diheds);

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
    typ1 = c->atoms[atm1].id;
    typ2 = c->atoms[atm2].id;
    typ3 = c->atoms[atm3].id;
    typ4 = c->atoms[atm4].id;
    fprintf(top, "%5d %5d %5d %5d 2             ; %s-%s-%s-%s\n", atm1, atm2,
        atm3, atm4, defs[typ4], defs[typ1], defs[typ2], defs[typ3]);
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
