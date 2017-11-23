#include "fragments.h"
#include "lammpstrj.h"
#include "xyz.h"

/* Reads a trajectory, checks if there are any CO or CO2 groups, and prints out
 * the id of C and O atoms if there are. 
 */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<1.8)||(typ>24&&dist<2.3))?1:0;
}


int main(int argc, char *argv[])
{
  int xyz = strcmp(strrchr(argv[1], '.'), ".xyz") != -1;

  FILE *f = fopen(argv[1], "r");

  int t, nat,
      elem[4] = {0,0,0,0}, // (Z-1)%4: H -> 0, C -> 1, N -> 2, O -> 3
      pos=0;

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

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    pos = 0;
    // read a frame
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
    // Do the work here!

    printf("10000 ");
    for (int i = 0; i < frg->nfrags; ++i)
    {
      memset(elem, 0, 4 * sizeof(int));
      for (int j = 0; j < frg->lfrags[i]; ++j)
      {
        elem[(c->atoms[frg->frags[pos]].Z-1)%4]++;
        pos++;
      }
      /* if(elem[1]==1&&(elem[3]==1||elem[3]==2)&&elem[0]==0&&elem[2]==0) */
      if(elem[1]==1&&elem[3]==2&&elem[0]==0&&elem[2]==0)
      {
        for (int j = pos-elem[3]-1; j < pos; ++j)
        {
          printf(" %d", frg->frags[j]);
        }
        printf("\n");
      }
    }
    /* printf("\n"); */


    // read the header information for the next frame
    if (xyz)
      XYZReadHeader(f, &nat, dm);
    else
      LMPReadHeader(f, &t, &nat, dm);

    // update cell size
    CrystalSetCell(c, dm);
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
