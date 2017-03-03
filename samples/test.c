#include "xyz.h"
#include "bonding.h"
#include <time.h>

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}

int oxidized(Crystal *c, BondingInfo *bnd, int index)
{
  for (int i = 0; i < bnd->nbonds[index]; ++i)
  {
    if (c->atoms[bnd->bonds[index][i]].Z == 8)
    {
      return 1;
    }
  }
  return 0;
}

int firstC(Crystal *c, BondingInfo *bnd, int index)
{
  for (int i = 0; i < bnd->nbonds[index]; ++i)
  {
    if (c->atoms[bnd->bonds[index][i]].Z == 6)
    {
      return bnd->bonds[index][i];
    }
  }
  return -1;
}

int main()
{
  FILE *f = fopen("0.xyz", "r");
  int nat, r, C, dir;
  double coor[3], dm[3];
  srand(time(NULL));

  XYZReadHeader(f, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  XYZReadFrame(f, c);
  CoarseBox *b = BoxInit(c, 2.0);
  BoxFill(c, b);
  BondingInfo *bnd = BondingInit(c);
  BondingPopulate(c, b, bnd, checkBonding);

  for (int i = 432; i < c->nat; ++i)
  {
    if ( c->atoms[i].Z == 6 && !oxidized(c, bnd, i) )
    {
      dir = pow(-1, i%2);
      r = rand() % 12;
      if (r == 1)
      {
        C = firstC(c, bnd, i);
        for (int j = 0; j < 3; ++j)
        {
          coor[j] = (c->atoms[i].coor[j] + c->atoms[C].coor[j])/2;
        }
        coor[2] += 1.3 * dir;
        if (b->binsn[BoxFromCoor(coor, b)] == 0)
        {
          printf("O %14.10f %14.10f %14.10f\n", coor[0], coor[1], coor[2]);
        }
      }
      else if (r == 2)
      {
        for (int j = 0; j < 3; ++j)
        {
          coor[j] = c->atoms[i].coor[j];
        }
        coor[2] += 1.8 * dir;
        if (b->binsn[BoxFromCoor(coor, b)] == 0)
        {
          printf("O %14.10f %14.10f %14.10f\n", coor[0], coor[1], coor[2]);
          printf("H %14.10f %14.10f %14.10f\n", coor[0], 0.7+coor[1], 0.2*dir+coor[2]);
        }
      }
    }
  }

  BoxDelete(b);
  BondingDelete(bnd);
  CrystalDelete(c);
  fclose(f);

  return 0;
}
