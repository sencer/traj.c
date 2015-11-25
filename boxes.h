#ifndef BONDING_H_5W4RONXR
#define BONDING_H_5W4RONXR

#include "crystal.h"
#include "util.h"

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
#endif /* end of include guard: BONDING_H_5W4RONXR */
