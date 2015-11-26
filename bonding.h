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
int BondingPopulate(Crystal *c, CoarseBox *box, BondingInfo *bnd);
int BondingClear(BondingInfo *bnd);
int checkBonding(Crystal *c, int i, int j);
int BondingPrint(BondingInfo *bnd);
int BondingFragments(BondingInfo *bnd);

#endif /* end of include guard: BONDING_H_SH2DFIOO */
