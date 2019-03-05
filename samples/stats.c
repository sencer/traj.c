// vi: fdm=syntax
#include "read.h"
#include "bonding.h"
#define PI 3.14159265
#define RAD2DEG 180/PI

void VecDiff(double u[3], double v[3], double diff[3])
{
  // difference of two 3D vectors
  for (int i = 0; i < 3; ++i)
  {
    diff[i] = u[i] - v[i];
  }
}

double VecLen2(double u[3])
{
  // length square of a 3D vector
  return pow(u[0], 2) + pow(u[1], 2) + pow(u[2], 2);
}

double Dist2(double u[3], double v[3])
{
  // squared euclidian distance between two 3D vectors
  double diff[3];
  VecDiff(u, v, diff);
  return VecLen2(diff);
}

double VecDot(double u[3], double v[3])
{
  // dot product of two 3D vectors
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

double VecAngle(double u[3], double v[3])
{
  // angle between two 3D vectors, in degrees
  double lu = sqrt(VecLen2(u)),
         lv = sqrt(VecLen2(v));
  return RAD2DEG * acos(VecDot(u, v) / (lu * lv));
}

int check_bonding(double len, int t1, int t2)
{
  int t = t1 + t2;
  return (len < 1.2 ||
      (len < 1.4 && t > 2) ||
      (len < 1.7 && t > 10) ||
      (len < 2.4 && t > 22));
}

int main(int argc, char *argv[])
{
  int t, xyz;
  FILE *f = open_file(argv[1], &xyz);
  Crystal *c = read_file(f, xyz, &t);
  CoarseBox *box = Box(c, 2.4);
  BondingInfo *bnd = Bonds(c, box, check_bonding, 0);

  int ia, ib, id, iframe = 0, eof = 0, typ;
  double r;
  Atom a, b, d;
  double u[3], v[3];

  FILE *bonds = fopen("bonds.dat", "w"),
       *angles = fopen("angles.dat", "w");

  do
  {
    for (int ia = 0; ia < c->nat; ++ia)
    {
      a = c->atoms[ia];
      if (a.Z == 1)
      {
        continue;
      }
      r = sqrt(pow(a.coor[0]-50, 2) + pow(a.coor[1]-50, 2) + pow(a.coor[2]-50, 2));
      typ = (a.Z == 22);
      for (int j = 0; j < bnd->nbonds[ia]; ++j)
      {
        ib = bnd->bonds[ia][j];
        b = c->atoms[ib];
        if ( (typ && b.Z != 8) || (!typ && b.Z != 22)) continue;
        VecDiff(b.coor, a.coor, u);
        fprintf(bonds, "%d %5.2f %5.3f\n", iframe, r, sqrt(VecLen2(u)));
        for (int k = j+1; k < bnd->nbonds[ia]; ++k)
        {
          id = bnd->bonds[ia][k];
          d = c->atoms[id];
          if ( (typ && d.Z != 8) || (!typ && d.Z != 22)) continue;
          VecDiff(d.coor, a.coor, v);
          fprintf(angles, "%d %d %5.2f %6.2f\n", iframe, typ, r, VecAngle(u, v));
        }
      }
    }
    iframe++;
    read_next(f, c, xyz, &t);

    BondingClear(bnd);
    BoxUpdate(c, box);
    BoxFill(c, box);
    BondingPopulate(c, box, bnd, check_bonding, 1);
    eof = fgetc(f);
    /* fseek(f, -1, SEEK_CUR); */
  }
  while(eof != EOF);

  fclose(f);
  fclose(bonds);
  fclose(angles);
  CrystalDelete(c);
  BondingDelete(bnd);
  BoxDelete(box);
  return 0;
}
