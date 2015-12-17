#include "unwrap.h"

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}


int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin,
       *f2= fopen("unwrap.lammpstrj", "w");
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
    BondingPopulate(c, box, bnd, checkBonding);
    // populate the fragments list
    BondingFragments(bnd);
    // Do the work here!
    BondingMergeFragments(c, bnd);
    LMPWriteFrame(f2, c, t);

    // read the header information for the next frame
    LMPReadHeader(f, &t, &nat, dm);
    // update cell size
    CrystalSetCell(c, dm);
    // Before moving to the next frame clear the bonding information
    BondingClear(bnd);
    // update the coarse graining box, in case cell dimensions changed
    BoxUpdate(c, box);
  }

  // free the memory
  fclose(f);
  fclose(f2);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
