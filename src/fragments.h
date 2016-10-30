#ifndef FRAGMENTS_H_DMNBQ35V
#define FRAGMENTS_H_DMNBQ35V
#include "bonding.h"

typedef struct fragment_info {
  int nat,            // number of atoms
      *frags,         // atoms listed in the order of fragments they belong
      *lfrags,        // length of each fragment
      nfrags;         // number of total fragments
} Fragments;

Fragments *FragmentsInit(BondingInfo *bnd);
void FragmentsDelete(Fragments *frg);
int FragmentsClear(Fragments *frg);
int FragmentsPopulate(BondingInfo *bnd, Fragments *frg);
int FragmentsMerge(Crystal *c, BondingInfo *bnd, Fragments *frg);

#endif /* end of include guard: FRAGMENTS_H_DMNBQ35V */
