#ifndef SURFACE_H_TPYQL5F6
#define SURFACE_H_TPYQL5F6

#include "vector.h"
#include <string.h>

typedef struct surface {
  double p[3], n[3];
} Surface;

void SurfaceCreate(int miller[3], double cell[3][3], double p[3], Surface *s);
double SurfaceDist(Surface *s, double p[3]);
void miller_to_plane(int miller[3], int v[3], int u[3]);
void miller_to_normal(int miller[3], double cell[3][3], double normal[3]);


#endif /* end of include guard: SURFACE_H_TPYQL5F6 */
