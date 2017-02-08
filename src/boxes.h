#ifndef BONDING_H_5W4RONXR
#define BONDING_H_5W4RONXR

#include "crystal.h"
#include "util.h"

#define MAPB 5       // Max Atoms Per Box TODO I should increase this, perhaps?
#define MBPA 12      // Max Bonds Per Atom

typedef struct coarse_box {
  int ntot,
      n2d,
      ngrid[3],
      (*bins)[MAPB],
      *binsn;
  double width, w[3];
} CoarseBox;

CoarseBox* BoxInit(Crystal *c, double width);
void BoxDelete(CoarseBox *box);
void BoxGetComponents(CoarseBox *box, int ind, int comp[3]);
int BoxGetIndice(CoarseBox *box, int comp[3]);
int BoxUpdate(Crystal *c, CoarseBox *box);
int BoxFill(Crystal *c, CoarseBox *box);
int BoxClear(CoarseBox *box);
int BoxFromCoor(double coor[3], CoarseBox *box);
void GetNeighboringBoxes(CoarseBox *box, int id, int ids[26]);
void BoxesOfAtoms(Crystal *c, CoarseBox *box, int *box_of_atom);
#endif /* end of include guard: BONDING_H_5W4RONXR */
