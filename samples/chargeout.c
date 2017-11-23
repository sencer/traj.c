#include "xyz.h"
#include "fragments.h"

/*
 * Calculates some specific species in a trajectory, like biphenyl groups
 * (defined very specifically for this system), CO2, O2 etc.
 */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0)||(typ>20&&dist<2.3))?1:0;
}


int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int nat, pos, types[3], tot_charge, o_charge;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  XYZReadHeader(f, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container
  FragmentsInfo *frg = FragmentsInit(bnd);

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    // read a frame
    XYZReadFrame(f, c);
    // assign the atoms to boxes
    BoxFill(c, box);
    // populate the bonding list
    BondingPopulate(c, box, bnd, checkBonding, 1);
    // populate the fragments list
    FragmentsPopulate(bnd, frg);
    // Do the printing here!

    pos = 0;
    printf("# new frame\n");
    for (int i = 0; i < frg->nfrags; ++i) {
      types[0] = 0;
      types[1] = 0;
      types[2] = 0;
      tot_charge = 0;
      for (int j = 0; j < frg->lfrags[i]; ++j) {
        tot_charge += c->atoms[frg->frags[pos]].id;
        if (c->atoms[frg->frags[pos]].Z == 1)
        {
          types[0]++;
        }
        else if (c->atoms[frg->frags[pos]].Z == 8)
        {
          types[1]++;
          o_charge = c->atoms[frg->frags[pos]].id;
        }
        else if (c->atoms[frg->frags[pos]].Z == 22)
        {
          types[2]++;
        }
        else { printf("%d has Z %d\n", frg->frags[pos], c->atoms[frg->frags[pos]].Z); }
        pos++;
      }
      if (types[0] == 2 && types[1] == 1 && types[2] == 0) // water
      {
        printf("%8.3f %8.3f\n", 1.0*o_charge/10000.0, 1.0*tot_charge/10000.0);
      }
    }
    printf("\n\n");
    fflush(stdout);

    // read the header information for the next frame
    XYZReadHeader(f, &nat, dm);
    // Before moving to the next frame clear the bonding information
    BondingClear(bnd);
    // and the fragmentation information
    FragmentsClear(frg);
    // update the coarse graining box, in case cell dimensions changed
    BoxUpdate(c, box);
  }

  // free the memory
  fclose(f);
  FragmentsDelete(frg);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
