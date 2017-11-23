#include "bonding.h"
#include "lammpstrj.h"
#include "xyz.h"

// Print out non-bonded H atoms

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (typ==2 && dist<1.0);
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int nat, t;
  double dm[3];

  int xyz = strcmp(strrchr(argv[1], '.'), ".xyz") != -1;
  // read the number of atoms and cell dimensions from lammpstrj
  if (xyz)
    XYZReadHeader(f, &nat, dm);
  else
    LMPReadHeader(f, &t, &nat, dm);

  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);

  // read a frame
  if (xyz)
    XYZReadFrame(f, c);
  else
    LMPReadFrame(f, c);
  // assign the atoms to boxes
  BoxFill(c, box);
  // populate the bonding list
  BondingPopulate(c, box, bnd, checkBonding, 1);

  Atom atm; // dbg
  for (int i = 0; i < c->nat; ++i)
  {
    /* if (c->atoms[i].Z == 1 && bnd->nbonds[i] == 0) */
    if (!(c->atoms[i].Z == 1 && bnd->nbonds[i] == 0)) //dbg
    {
      /* printf("%s:%d\n", argv[1], i+3); */
      atm = c->atoms[i]; //dbg
      printf("%s %10.6f %10.6f %10.6f\n", PT_Symbol(atm.Z), atm.coor[0], atm.coor[1], atm.coor[2]); //dbg
    }
  }

  // free the memory
  fclose(f);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
