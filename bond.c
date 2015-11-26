#include "bond.h"

int main(int argc, char *argv[])
{
  FILE *f = fopen(argv[1], "r");
  int t, nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    // read a frame
    LMPReadFrame(f, c);
    // assign the atoms to boxes
    BoxFill(c, box);
    // populate the bonding list
    BondingPopulate(c, box, bnd);
    // populate the fragments list
    BondingFragments(bnd);
    // Do the printing here!

    // read the header information for the next frame
    LMPReadHeader(f, &t, &nat, dm);
    // Before moving to the next frame clear the bonding information
    BondingClear(bnd);
    // update the coarse graining box, in case cell dimensions changed
    BoxUpdate(c, box);
  }

  // free the memory
  fclose(f);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
