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

  // we will use up to four atoms (when writing Dihedrals for example)
  // so lets allocate some memory for them.
  int typ1, typ2, typ3, typ4;

  char
    // we need to teach opls atom types for each of these atoms
    opls_type[10][9] = {
      "opls_147", "opls_166", "opls_182", "opls_267", "opls_167", "opls_179",
      "opls_269", "opls_268", "opls_168", "opls_270"
    },
    // we will also want to output the atom types in .top file for human reader
    defs[10][4] = {
      "CA", "CF", "CT", "C", "OH", "OS", "O", "OH2", "HO", "HO2"
    };

  double
    // charges of each atom type, taken from literature TODO where exactly?
    charges[10] = {
      0.0,  0.15,  0.14,  0.52, -0.585, -0.28, -0.44, -0.53, 0.435, 0.45
    },
    // let's keep the total charge to make sure the whole molecule is neutral
    tot_charge,
    // masses of each atom type
    masses[10] = {
      12.01100, 12.01100, 12.01100, 12.01100, 15.99940, 15.99940, 15.99940,
      15.99940, 1.00794,  1.00794
    };

  // Read the xyz file to a Crystal
  Crystal *c = ReadXYZ(argv[1]);
  ConvertAtomIDsToTypes(c);

  // Open the argv[2] to write the .top file
  FILE *top = fopen(argv[2], "w");

  // Now let's write the [atoms] section
  fprintf(top, "\n[ atoms ]\n");
  for (int i = 0; i < c->nat; ++i)
  {
    typ1 = c->atoms[i].id;
    fprintf(top, "%5d    %-10s   1     GO_s      %s %4d %9.5f %9.5f\n",
        i+1, opls_type[typ1], PT_Symbol(c->atoms[i].Z), i+1, charges[typ1],
        masses[typ1]);
  }

  return 0;
}

// vi: fdm=syntax
