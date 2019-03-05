#ifndef VECTOR_H_GKHNJDFW
#define VECTOR_H_GKHNJDFW
#include <math.h>

double VecDot(double v[3], double u[3]);
void VecCross(double v[3], double u[3], double ret[3]);
void VecNormalize(double v[3]);
void VecInvert(double v[3]);
void MatDot(double mat[3][3], int v[3], double ret[3]);

int gcd(int a, int b);
int lcm(int a, int b);

#endif /* end of include guard: VECTOR_H_GKHNJDFW */
