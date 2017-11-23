#include "xyz.h"
#include "topo.h"

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}

int main(int argc, char *argv[])
{
  FILE *f = fopen(argv[1], "r");
  int nat;
  double dm[3];

  XYZReadHeader(f, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  XYZReadFrame(f, c);
  fclose(f);

  CoarseBox *box = BoxInit(c, 2.0);
  BoxFill(c, box);

  BondingInfo *bnd = BondingInit(c);
  BondingPopulate(c, box, bnd, checkBonding, 1);
  BoxDelete(box);

  int nangles;
  int *angles = Angles(c, bnd, &nangles);

  for (int i = 0; i < nangles; ++i)
  {
    printf("%d %d %d\n", angles[3*i], angles[3*i+1], angles[3*i+2]);
  }
  free(angles);

  int ndiheds;
  int *diheds = Dihedrals(c, bnd, &ndiheds);
  for (int i = 0; i < ndiheds; ++i)
  {
    printf("%d %d %d %d\n", diheds[4*i], diheds[4*i+1], diheds[4*i+2], diheds[4*i+3]);
  }
  free(diheds);

  int nimps;
  int *imps = Impropers(c, bnd, &nimps);
  for (int i = 0; i < nimps; ++i)
  {
    printf("%d %d %d %d\n", imps[4*i], imps[4*i+1], imps[4*i+2], imps[4*i+3]);
  }
  free(imps);

  BondingDelete(bnd);
  CrystalDelete(c);

  return 0;
}
