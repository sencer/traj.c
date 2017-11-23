#include "lammpstrj.h"
#include "fragments.h"

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>23 && dist<2.4))?1:0;
}

int other(Crystal *c, FragmentsInfo *frg)
{
  int form[3], fpos = 0, h=0, o=0, oh=0, h2=0, o2=0, h2o=0, h2o2=0;
  double tio[3] = {0, 0, 0}; 

  for (int i = 0; i < frg->nfrags; ++i)
  {
    memset(form, 0, 3*sizeof(int));
    for (int j = 0; j < frg->lfrags[i]; ++j)
    {
      switch (c->atoms[frg->frags[fpos++]].Z) {
        case 1:
          form[0]++;
          break;
        case 8:
          form[1]++;
          break;
        case 22:
          form[2]++;
          break;
      }
    }
    if (form[2] == 0)
    {
      if      (form[0] == 1 && form[1] == 0) { h++; }
      else if (form[0] == 2 && form[1] == 0) { h2++; }
      else if (form[0] == 2 && form[1] == 1) { h2o++; }
      else if (form[0] == 1 && form[1] == 1) { oh++; }
      else if (form[0] == 0 && form[1] == 1) { o++; }
      else if (form[0] == 0 && form[1] == 2) { o2++; }
      else if (form[0] == 2 && form[1] == 2) { h2o2++; }
    }
    else
    {
      tio[0] += form[2];
      tio[1] += form[1];
      tio[2] += form[0];
    }
  }
  printf("TiO %6.4f H %6.4f H %d H2 %d O %d OH %d O2 %d H2O %d H2O2 %d Vol %12.6f\n", tio[1]/tio[0],
      tio[2]/tio[0], h, h2, o, oh, o2, h2o, h2o2,c->dm[0]*c->dm[1]);
  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.4);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container
  FragmentsInfo *frg = FragmentsInit(bnd);

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    // read a frame
    LMPReadFrame(f, c);
    // assign the atoms to boxes
    BoxFill(c, box);
    // populate the bonding list
    BondingPopulate(c, box, bnd, checkBonding, 1);
    // populate the fragments list
    FragmentsPopulate(bnd, frg);
    // Do the printing here!
    printf("%-8d  ", t);
    other(c, frg);
    fflush(stdout);

    // read the header information for the next frame
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
