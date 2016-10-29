#ifndef CRYSTAL_H_F4OAZMJQ
#define CRYSTAL_H_F4OAZMJQ
#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef struct atom {
  int Z;                           // atomic number
  double coor[3];                  // coordinates
  int id;
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
#endif /* end of include guard: CRYSTAL_H_F4OAZMJQ */
