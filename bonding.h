#ifndef BONDING_H_SH2DFIOO
#define BONDING_H_SH2DFIOO
#include <stdio.h>
#include "crystal.h"
#include "boxes.h"


typedef struct bonding_info {
  int nat, 
      (*bonds)[MBPA],
      *bondsn,
      *lfrags,
      *frags,
      nfrags;
} BondingInfo;

BondingInfo *BondingInit(Crystal *c);
void BondingDelete(BondingInfo *bnd);
int BondingPopulate(Crystal *c, CoarseBox *box, BondingInfo *bnd, int (*cb)(double, int, int));
int BondingDoBonding(Crystal *c, BondingInfo *bnd, int atm1, int atm2, int (*cb)(double, int, int));
int BondingClear(BondingInfo *bnd);
int BondingPrint(BondingInfo *bnd);
int BondingFragments(BondingInfo *bnd);
int BondingMergeFragments(Crystal *c, BondingInfo *bnd);

#endif /* end of include guard: BONDING_H_SH2DFIOO */
