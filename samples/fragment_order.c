#include "xyz.h"
#include "fragments.h"

/* Orders the atoms in a frame so that atoms in each fragment are sequential */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<1.8)||(typ>24&&dist<2.3))?1:0;
}


int main(int argc, char *argv[])
{
  FILE *f = fopen(argv[1], "r");
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
  // ...and a fragmentation information container
  FragmentsInfo *frg = FragmentsInit(bnd);

  XYZReadFrame(f, c);
  BoxFill(c, box);
  BondingPopulate(c, box, bnd, checkBonding, 1);
  FragmentsPopulate(bnd, frg);

  int cur;
  for (int i = 0; i < c->nat; ++i)
  {
    cur = frg->frags[i];
    printf("%-2s %10.6f %10.6f %10.6f\n", PT_Symbol(c->atoms[cur].Z),
      c->atoms[cur].coor[0], c->atoms[cur].coor[1], c->atoms[cur].coor[2]);
  }

  // free the memory
  fclose(f);
  FragmentsDelete(frg);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
