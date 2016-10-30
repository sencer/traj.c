#ifndef BONDING_H_SH2DFIOO
#define BONDING_H_SH2DFIOO
#include <stdio.h>
#include <math.h>
#include "crystal.h"
#include "boxes.h"

// Defines a new type, BondingInfo, that holds a table of bonding for
// a Crystal, and information on the fragmentation structure - by which I mean,
// the information of which atom belongs to which molecule/fragment.

// Also provides methods to fill these information, and use this information
// to merge fragments broken because of periodic boundary conditions.

typedef struct bonding_info {
  int nat,            // number of atoms
      (*bonds)[MBPA], // bonding table, for each atom max MBPA other listed
      *bondsn;        // number of bonds of each atom, since often it's not MBPA
} BondingInfo;

// Allocates the memory needed for bonding info
BondingInfo *BondingInit(Crystal *c);
// De-allocates the memory used for bonding info
void BondingDelete(BondingInfo *bnd);
// Fills in the bonding information. Needs:
// a Crystal
// a CoarseBox -see the header file for its explanation-
// a BondingInfo to fill the data in (previously initialized), and
// a Function that receives a distance, and atomic numbers of 2 atoms
// This function should return 0 for not bonded, and 1 for bonded.
int BondingPopulate(Crystal *c, CoarseBox *box, BondingInfo *bnd, int (*cb)(double, int, int));
// Clear the bonding information, keep the memory. This needs to be used
// Because part of the table will not be always filled.
int BondingClear(BondingInfo *bnd);

int BondingDoBonding(Crystal *c, BondingInfo *bnd, int atm1, int atm2, int (*cb)(double, int, int));
int BondingPrint(BondingInfo *bnd);

#endif /* end of include guard: BONDING_H_SH2DFIOO */
