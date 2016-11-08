#include "bonding.h"
#include "xyz.h"

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (typ==2 &&dist<1.2)?1:0;
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  XYZReadHeader(f, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);

  // read a frame
  XYZReadFrame(f, c);
  // assign the atoms to boxes
  BoxFill(c, box);
  // populate the bonding list
  BondingPopulate(c, box, bnd, checkBonding);

  for (int i = 0; i < c->nat; ++i)
  {
    if (c->atoms[i].Z == 1 && bnd->bondsn[i] == 0)
    {
      printf("%s:%d\n", argv[1], i+3);
    }
  }

  // free the memory
  fclose(f);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
