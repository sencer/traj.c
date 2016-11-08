#include "lammpstrj.h"
#include "crystal.h"
#define LIM 20.0 // An upper-limit for the longest possible travel per step.

int main(int argc, char *argv[])
{
  FILE *f = fopen(argv[1], "r");
  FILE *o = fopen("wrapped.lammpstrj", "w");
  int t, nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);

  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    LMPReadFrame(f, c);
    LMPWriteFrame(o, c, t);
    LMPReadHeader(f, &t, &nat, dm);
    CrystalSetCell(c, dm);
  }

  // free the memory
  fclose(f);
  fclose(o);
  CrystalDelete(c);

  return 0;
}
