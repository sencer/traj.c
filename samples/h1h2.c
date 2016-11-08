#include "bonding.h"
#include "xyz.h"
#define START 7539
#define NUM 762

int checkBonding(double dist, int t1, int t2)
{
  return (t1+t2==2&&dist<1.4)?1:0;
}

int main(int argc, char *argv[])
{
  int h_start = START, h_num = NUM;
  if (argc < 2)
  {
    printf("Please give a file name.\n");
    return -1;
  }

  FILE *f = fopen(argv[1], "r");
  int nat, *list;
  double dm[3];
  XYZReadHeader(f, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  CoarseBox *box = BoxInit(c, 2.0);
  XYZReadFrame(f, c);
  fclose(f);
  BondingInfo *bnd = BondingInit(c);
  BoxFill(c, box);
  BondingPopulate(c, box, bnd, checkBonding);

  list = malloc(h_num*sizeof(int));
  memset(list, 0, h_num*sizeof(int));

  for (int i = h_start; i < h_start+h_num; ++i)
  {
    if (bnd->bondsn[i] > 1)
    {
      printf("Hydrogen %d (line: %d) has %d bonds!\n", i-h_start+1, i+3, bnd->bondsn[i]);
    }
    else if(bnd->bondsn[i] == 0)
    {
      printf("Hydrogen %d is not bonded!\n", i);
      fprintf(stderr, "%s:%d\n", argv[1], 3+i);
    }
    else
    {
      // TODO this uses a temprorary hack in the xyz.c, which translates 
      // Z = 0 to H1, and Z = 2 to H2 rather than He!
      if (!list[i-h_start])
      {
        c->atoms[i].Z = 0;
      }
      if (!list[bnd->bonds[i][0]-h_start])
      {
        c->atoms[bnd->bonds[i][0]].Z = 2;
      }
      list[i-h_start] = 1;
      list[bnd->bonds[i][0]-h_start] = 1;
    }
  }
  f = fopen("out.xyz", "w");
  free(list);
  XYZWriteFrame(f, c);
  fclose(f);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);
  return 0;
}
