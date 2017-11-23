#ifndef FRAGMENTS_H_DMNBQ35V
#define FRAGMENTS_H_DMNBQ35V
#include "bonding.h"

typedef struct fragment_info {
  int nat,            // number of atoms
      *frags,         // atoms listed in the order of fragments they belong
      *lfrags,        // length of each fragment
      nfrags;         // number of total fragments
} FragmentsInfo;

FragmentsInfo *FragmentsInit(BondingInfo *bnd);
void FragmentsDelete(FragmentsInfo *frg);
int FragmentsClear(FragmentsInfo *frg);
int FragmentsPopulate(BondingInfo *bnd, FragmentsInfo *frg);
int FragmentsMerge(Crystal *c, BondingInfo *bnd, FragmentsInfo *frg);
int* FragmentHead(FragmentsInfo *frg, int n);
int* FragmentList(FragmentsInfo *frg);
FragmentsInfo *Fragments(BondingInfo *bnd);

#endif /* end of include guard: FRAGMENTS_H_DMNBQ35V */
