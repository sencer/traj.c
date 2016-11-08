#include "lammpstrj.h"
#include "crystal.h"
#define LIM 20.0 // An upper-limit for the longest possible travel per step.

void SwapAtoms(Crystal *c1, Crystal *c2)
{
  void* tmp = c1->atoms;
  c1->atoms = c2->atoms;
  c2->atoms = tmp;
}

int main(int argc, char *argv[])
{
  FILE *f = fopen(argv[1], "r");
  FILE *o = fopen("travel.lammpstrj", "w");
  int t, nat;
  double dm[3], (*shift)[3];

  LMPReadHeader(f, &t, &nat, dm);

  shift  = malloc(nat * sizeof(double[3]));
  memset(shift, 0, nat*3*sizeof(double));

  Crystal *prev = CrystalInit(nat, dm);
  LMPReadFrame(f, prev);

  Crystal *cur = CrystalInit(nat, dm);
  Crystal *out = CrystalInit(nat, dm);
  memcpy(out->atoms, prev->atoms, nat*sizeof(Atom));

  while(!feof(f))
  {
    LMPWriteFrame(o, out, t);
    LMPReadHeader(f, &t, &nat, dm);

    CrystalSetCell(cur, dm);
    CrystalSetCell(out, dm);
    LMPReadFrame(f, cur);


    for (int i = 0; i < nat; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        shift[i][j] += round((prev->atoms[i].coor[j]-cur->atoms[i].coor[j])/dm[j])*dm[j];
        out->atoms[i].coor[j] = cur->atoms[i].coor[j] + shift[i][j];
      }
    }
    SwapAtoms(prev, cur);
  }


  free(shift);
  fclose(f);
  fclose(o);
  CrystalDelete(prev);
  CrystalDelete(cur);
  CrystalDelete(out);

  return 0;
}
