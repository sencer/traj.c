// vi: fdm=syntax
#include "lammpstrj.h"
#include "fragments.h"

/*
 * Calculates some specific species in a trajectory, like biphenyl groups
 * (defined very specifically for this system), CO2, O2 etc.
 */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}

void biphenyl(int nat, FragmentsInfo *frg)
{
  int l[12] = {3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15},
      ldiff[nat], lend,
      pfrag,
      cnt = nat/38, icnt = cnt;

  for (int i = 0; i < icnt; ++i)
  {
    // compare `l` with all fragments. an awful implementation.
    // TODO it shouldn't check a fragment for more than one biphenyls
    pfrag = 0;
    for (int j = 0; j < frg->nfrags; ++j)
    {
      listDiff(l, 12, frg->frags+pfrag, frg->lfrags[j], ldiff, &lend);
      if(lend == 0) // this means all the atoms are in a single fragment
      {
        break;
      }
      else if (lend < 12) // this means atoms are fragmented!
      {
        cnt--;
        break;
      }
      pfrag += frg->lfrags[j];
    }
    for (int j = 0; j < 12; ++j) { l[j] += 38; }
  }
  printf("%4d", cnt);
  /* return cnt; */
}

int other(Crystal *c, FragmentsInfo *frg)
{
  int form[4], fpos = 0, co = 0, co2 = 0, n2 = 0, n = 0, h2 = 0, o2 = 0, h2o=0, oh=0;

  for (int i = 0; i < frg->nfrags; ++i)
  {
    memset(form, 0, 4*sizeof(int));
    for (int j = 0; j < frg->lfrags[i]; ++j)
    {
      switch (c->atoms[frg->frags[fpos++]].Z) {
        case 1:
          form[0]++;
          break;
        case 6:
          form[1]++;
          break;
        case 7:
          form[2]++;
          break;
        case 8:
          form[3]++;
          break;
      }
    }
    if      (form[2] == 1 && form[0]+form[1]+form[3]==0) { n++; }
    else if (form[2] == 2 && form[0]+form[1]+form[3]==0) { n2++; }
    else if (form[1] == 1 && form[3] == 1 && form[0]+form[2]==0) { co++; }
    else if (form[1] == 1 && form[3] == 2 && form[0]+form[2]==0) { co2++; }
    else if (form[0] == 2 && form[1]+form[2]+form[3]==0) { h2++; }
    else if (form[0] == 2 && form[3] == 1 && form[1]+form[2]==0) { h2o++; }
    else if (form[0] == 1 && form[3] == 1 && form[1]+form[2]==0) { oh++; }
    else if (form[3] == 2 && form[0]+form[1]+form[2]==0) { o2++; }
  }
  printf("%5d%5d%5d%5d%5d%5d%5d%5d%5d\n", n2, n, co, co2, h2, o2, h2o, oh, frg->nfrags);
  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat;
  double dm[3];
  printf("# Time      BP   N2    N   CO  CO2   H2   O2   H2O  OH  NF\n");

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);
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
    LMPReadFrame(f, c);
    // assign the atoms to boxes
    BoxFill(c, box);
    // populate the bonding list
    BondingPopulate(c, box, bnd, checkBonding, 1);
    // populate the fragments list
    FragmentsPopulate(bnd, frg);
    // Do the printing here!
    printf("%-8d  ", t);
    biphenyl(bnd->nat, frg);
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
