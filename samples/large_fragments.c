#include "fragments.h"
#include "lammpstrj.h"
#include "xyz.h"
#define LARGE_FRAG_SIZE 100

/* Print out only the fragments larger than LARGE_FRAG_SIZE atoms as separate
 * xyz files */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<1.8)||(typ>24&&dist<2.3))?1:0;
}


int main(int argc, char *argv[])
{
  int xyz = strcmp(strrchr(argv[1], '.'), ".xyz") != -1;

  FILE *f = fopen(argv[1], "r");

  int t, nat, pos=0;

  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  if (xyz)
    XYZReadHeader(f, &nat, dm);
  else
    LMPReadHeader(f, &t, &nat, dm);

  // read the number of atoms and cell dimensions from lammpstrj
  /* LMPReadHeader(f, &t, &nat, dm); */
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container
  FragmentsInfo *frg = FragmentsInit(bnd);

  if (xyz)
    XYZReadFrame(f, c);
  else
    LMPReadFrame(f, c);

  // assign the atoms to boxes
  BoxFill(c, box);
  // populate the bonding list
  BondingPopulate(c, box, bnd, checkBonding, 1);
  // populate the fragments list
  FragmentsPopulate(bnd, frg);

  // read the header information for the next frame
  if (xyz)
    XYZReadHeader(f, &nat, dm);
  else
    LMPReadHeader(f, &t, &nat, dm);

  FILE *fragment;
  char fname[20];
  int j = 0;
  Atom atm;
  for (int i = 0; i < frg->nfrags; ++i)
  {
    if(frg->lfrags[i] > LARGE_FRAG_SIZE)
    {
      sprintf(fname, "fragment_%d.xyz", j++);
      fragment = fopen(fname, "w");
      fprintf(fragment, "%d\ncelldm %f %f %f\n", frg->lfrags[i], c->dm[0], c->dm[1], c->dm[2]);
      for (int k = 0; k < frg->lfrags[i]; ++k)
      {
        atm = c->atoms[frg->frags[pos++]];
        if (atm.id > 0)
        {
          fprintf(fragment, "%2s%-3d", PT_Symbol(atm.Z), atm.id);
        }
        else
        {
          fprintf(fragment, " %-4s", PT_Symbol(atm.Z));
        }
        fprintf(fragment, " %12.9f %12.9f %12.9f\n", atm.coor[0],
            atm.coor[1],
            atm.coor[2]);
      }
      fclose(fragment);
    }
    else
    {
      pos += frg->lfrags[i];
    }
  }

  // free the memory
  fclose(f);
  FragmentsClear(frg);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
