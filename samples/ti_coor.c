// vi: fdm=syntax
#include "lammpstrj.h"
#include "xyz.h"
#include "bonding.h"

/*
 * Calculates some specific species in a trajectory, like biphenyl groups
 * (defined very specifically for this system), CO2, O2 etc.
 */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>23 && dist<2.6))?1:0;
}

int main(int argc, char *argv[])
{
  FILE *f = fopen(argv[1], "r");
  FILE *o = fopen(argv[2], "w");
  int t, nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.6);
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
    BondingPopulate(c, box, bnd, checkBonding, 1);

    // For each atom
    for (int i = 0; i < c->nat; ++i)
    {
      c->atoms[i].id = bnd->nbonds[i];
      /* if (c->atoms[i].Z == 8) */
      /* { */
      /*   c->atoms[i].id = 3-bnd->nbonds[i]; */
      /* } */
      /* else if (c->atoms[i].Z == 22) */
      /* { */
      /*   c->atoms[i].id = 6-bnd->nbonds[i]; */
      /* } */
      /* else{ */
      /*   c->atoms[i].id = 0; */
      /* } */
    }

    XYZWriteFrame(o, c);
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
  fclose(o);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
