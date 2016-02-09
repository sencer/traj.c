#include "bond.h"

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}

void biphenyl(int nat, Fragments *frg)
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

int other(Crystal *c, Fragments *frg)
{
  int form[4], fpos = 0, co = 0, n2 = 0, n = 0, h2 = 0;

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
    else if (form[0] == 2 && form[1]+form[2]+form[3]==0) { h2++; }
    else if (form[2] == 2 && form[0]+form[1]+form[3]==0) { n2++; }
    else if (form[1] == 1 && form[3] == 1 && form[0]+form[2]==0) { co++; }
  }
  printf("%5d%5d%5d%5d%5d\n", n2, n, co, h2, frg->nfrags);
  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat;
  double dm[3];
  printf("# Time      BP   N2    N   CO   H2   NF\n");

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container
  Fragments *frg = FragmentsInit(bnd);

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
