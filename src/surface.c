#include "surface.h"

void miller_to_plane(int miller[3], int v[3], int u[3])
{
  int izero[3], nzero = 0, tmp;

  // store the indices of zero elements from the beginning of izero array the
  // elements after nzero will be nonzero elements
  for (int i = 0; i < 3; ++i)
  {
    v[i] = 0;
    u[i] = 0;
    if (miller[i] == 0)
      izero[nzero++] = i;
    else
      izero[2-i+nzero] = i;
  }

  if (nzero == 2)
  {
    v[izero[0]] = 1;
    u[izero[1]] = 1;
  }
  else if (nzero == 1)
  {
    v[izero[0]] = 1;
    tmp = lcm(miller[izero[1]], miller[izero[2]]);
    u[izero[1]] = -tmp / miller[izero[1]];
    u[izero[2]] =  tmp / miller[izero[2]];
  }
  else
  {
    tmp = lcm(miller[izero[0]], miller[izero[1]]);
    v[0] = -tmp / miller[0];
    v[1] =  tmp / miller[1];
    tmp = lcm(miller[izero[0]], miller[izero[2]]);
    u[0] = -tmp / miller[0];
    u[2] =  tmp / miller[2];
  }
}

void miller_to_normal(int miller[3], double cell[3][3], double normal[3])
{
  int v[3], u[3];
  double V[3], U[3];
  miller_to_plane(miller, v, u);
  MatDot(cell, v, V);
  MatDot(cell, u, U);
  VecCross(U, V, normal);
  VecNormalize(normal);
  if( miller[0]*normal[0] + miller[1]*normal[1] + miller[2]*normal[2] < 0)
  {
    VecInvert(normal);
  }
}

void SurfaceCreate(int miller[3], double cell[3][3], double p[3], Surface *s)
{
  memcpy(s->p, p, 3*sizeof(double));
  miller_to_normal(miller, cell, s->n);
}

double SurfaceDist(Surface *s, double p[3])
{
  double dp[3] = { p[0]-s->p[0], p[1]-s->p[1], p[2]-s->p[2] };
  return VecDot(dp, s->n);
}
