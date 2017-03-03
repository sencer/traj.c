#ifndef TOPO_H_YBVO0NUX
#define TOPO_H_YBVO0NUX
#include "bonding.h"

int *Angles(Crystal *c, BondingInfo *bnd,  int *n);
int *Dihedrals(Crystal *c, BondingInfo *bnd, int *n);
int *Impropers(Crystal *c, BondingInfo *bnd,  int *n);

#endif /* end of include guard: TOPO_H_YBVO0NUX */
