#include "pdb.h"
#include "xyz.h"
#include "crystal.h"

/* convert pdb to lammpstrj */

int main(int argc, char *argv[])
{

  FILE *pdb = fopen(argv[1], "r"),
       *lmp = fopen(argv[2], "w");

  int t=0, nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  XYZReadHeader(pdb, &nat, dm);

  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);

  while(!feof(pdb))
  {
    XYZReadFrame(pdb, c);
    PDBWriteFrame(lmp, c);
    XYZReadHeader(pdb, &nat, dm);
    CrystalSetCell(c, dm);
  }

  // free the memory
  fclose(pdb);
  fclose(lmp);
  CrystalDelete(c);

  return 0;
}
