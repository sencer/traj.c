#include "crystal.h"

Crystal *CrystalInit(int nat, double dm[3])
{
  Crystal *c = malloc(sizeof(Crystal));
  CrystalSetNat(c, nat);
  CrystalSetCell(c, dm);
  return c;
}

void CrystalSetNat(Crystal *c, int nat)
{
  c->nat = nat;
  c->atoms = malloc(nat * sizeof(Atom));
  memset(c->atoms, 0, nat * sizeof(Atom));
}

void CrystalSetCell(Crystal *c, double dm[3])
{
  memcpy(c->dm, dm, 3 * sizeof(double));
}

void CrystalDelete(Crystal *c)
{
  free(c->atoms);
  free(c);
}

double CrystalDist(Crystal *c, double c1[3], double c2[3])
{
  double dist, tot = 0, *dm = c->dm;

  for (int i = 0; i < 3; ++i)
  {
    dist = 0.5 * (dm[i] - fabs(dm[i] - 2*fabs(c1[i] - c2[i])));
    tot += dist*dist;
  }

  return sqrt(tot);
}

