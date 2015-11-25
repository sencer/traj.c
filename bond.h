#ifndef BOND_H_CTJX7VU0
#define BOND_H_CTJX7VU0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int isInList(int elem, int *list, int llen);
int listDiff(int *l1, int len1, int *l2, int len2, int *ld, int *lend);
int mod(int x, int y);

typedef struct atom {
  int Z;                           // atomic number
  double coor[3];                  // coordinates
} Atom;

typedef struct crystal {
  int nat;                         // number of atoms
  double dm[3];                    // cell dimensions
  Atom *atoms;                     // atoms
} Crystal;

Crystal *CrystalInit(int nat, double dm[3]);
void CrystalSetNat(Crystal *c, int nat);
void CrystalSetCell(Crystal *c, double dm[3]);
void CrystalDelete(Crystal *c);
double CrystalDist(Crystal *c, int a, int b);

#define WIDTH 2.0    // Coarse Grained Box Size in Angstroms
#define MAPB 5       // Max Atoms Per Box TODO I should increase this, perhaps?
#define MBPA 12      // Max Bonds Per Atom

typedef struct coarse_box {
  int ntot,
      n2d,
      ngrid[3],
      (*bins)[MAPB],
      *binsn;

} CoarseBox;

CoarseBox* BoxInit(Crystal *c);
void BoxDelete(CoarseBox *box);
void BoxGetComponents(CoarseBox *box, int ind, int comp[3]);
int BoxGetIndice(CoarseBox *box, int comp[3]);
int BoxUpdate(Crystal *c, CoarseBox *box);
int BoxFill(Crystal *c, CoarseBox *box);
int BoxClear(CoarseBox *box);

typedef struct bonding_info {
  int nat, 
      (*bonds)[MBPA],
      *bondsn,
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

int LMPReadHeader(FILE *f, int *t, int *nat, double *dm);
int LMPReadFrame(FILE *f, Crystal *c);

#endif /* end of include guard: BOND_H_CTJX7VU0 */
