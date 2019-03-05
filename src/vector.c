#include "vector.h"

double VecDot(double v[3], double u[3])
{
  return v[0]*u[0] + v[1]*u[1] + v[2]*u[2];
}

void VecCross(double v[3], double u[3], double ret[3])
{
  ret[0] = v[1]*u[2] - v[2]*u[1];
  ret[1] = v[2]*u[0] - v[0]*u[2];
  ret[2] = v[0]*u[1] - v[1]*u[0];
}

void VecNormalize(double v[3])
{
  double len = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  for (int i = 0; i < 3; ++i)
  {
    v[i] /= len;
  }
}

void MatDot(double mat[3][3], int v[3], double ret[3])
{
  for (int i = 0; i < 3; ++i)
  {
    ret[i] = mat[i][0]*v[0] + mat[i][1]*v[1] + mat[i][2]*v[2];
  }
}

void VecInvert(double v[3])
{
  for (int i = 0; i < 3; ++i)
  {
    v[i] = -v[i];
  }
}

int gcd(int a, int b)
{
  if (a == 0)
  {
    return b;
  }
  return gcd(b % a, a);
}

int lcm(int a, int b)
{
  return (a*b)/gcd(a, b);
}
