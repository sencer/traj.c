#include "bonding.h"

BondingInfo *BondingInit(Crystal *c)
{
  BondingInfo *bnd = malloc(sizeof(BondingInfo));

  bnd->nat = c->nat;
  bnd->bonds  = malloc(bnd->nat * sizeof(int[MBPA]));
  bnd->bondsn = malloc(bnd->nat * sizeof(int));

  BondingClear(bnd);

  return bnd;
}

void BondingDelete(BondingInfo *bnd)
{
  free(bnd->bonds);
  free(bnd->bondsn);
  free(bnd);
}

int BondingDoBonding(Crystal *c, BondingInfo *bnd, int i, int j, int (*cb)(double, int, int))
{
  double dist = CrystalDist(c, i, j);
  if(cb(dist, c->atoms[i].Z, c->atoms[j].Z))
  {
    bnd->bonds[i][bnd->bondsn[i]++] = j;
    bnd->bonds[j][bnd->bondsn[j]++] = i;
  }
  return 0;
}

int BondingPopulate(Crystal *c, CoarseBox *box, BondingInfo *bnd, int (*cb)(double, int, int))
{
  int icomp[3], jcomp[3], jbox, atm1, atm2;

  // for each box
  for (int ibox = 0; ibox < box->ntot; ++ibox)
  {
    // if there are atoms in the box
    if (box->binsn[ibox])
    {
      // cross-check bonding between the atoms in box ibox
      for (int i = 0; i < box->binsn[ibox]-1; ++i)
      {
        atm1 = box->bins[ibox][i];
        for (int j = i+1; j < box->binsn[ibox]; ++j)
        {
          atm2 = box->bins[ibox][j];
          BondingDoBonding(c, bnd, atm1, atm2, cb);
        }
      }
      // now check the neighboring boxes - only the 13 of them in (+) direction
      BoxGetComponents(box, ibox, icomp);
      for (int j = 14; j < 27; ++j)
      {
        jcomp[0] = icomp[0] + j / 9 - 1;
        jcomp[1] = icomp[1] + (j % 9) / 3 - 1;
        jcomp[2] = icomp[2] + j % 3 - 1;
        jbox = BoxGetIndice(box, jcomp);
        // for each atom in box jbox check bonding
        for (int k = 0; k < box->binsn[jbox]; ++k)
        {
          atm1 = box->bins[jbox][k];
          for (int i = 0; i < box->binsn[ibox]; ++i)
          {
            atm2 = box->bins[ibox][i];
            BondingDoBonding(c, bnd, atm1, atm2, cb);
          }
        }
      }
    }
  }
  return 0;
}

int BondingClear(BondingInfo *bnd)
{
  memset(bnd->bondsn, 0, bnd->nat * sizeof(int));
  return 0;
}

int BondingPrint(BondingInfo *bnd)
{
  for (int i = 0; i < bnd->nat; ++i)
  {
    printf("%-4d: ", i);
    for (int j = 0; j < bnd->bondsn[i]; ++j)
    {
      printf("%d ", bnd->bonds[i][j]);
    }
    printf("\n");
  }
  return 0;
}