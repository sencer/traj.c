#include <stdlib.h>
#include "lammpstrj.h"

Crystal *ReadFrame(FILE *f)
{
  int t, nat;
  double dm[3];

  LMPReadHeader(f, &t, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);

  LMPReadFrame(f, c);

  return c;
}

int main(int argc, char *argv[])
{
  FILE *f0 = fopen(argv[1], "r"),
       *f1 = fopen(argv[2], "r");

  int nimages = atoi(argv[3]);
  double c = 1.0/(nimages+1), (*diff)[3];

  Crystal *first = ReadFrame(f0),
          *last  = ReadFrame(f1);

  char img[30];

  fclose(f0);
  fclose(f1);

  diff = malloc(first->nat*sizeof(double[3]));

  for (int i = 0; i < first->nat; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      diff[i][j] = c*(last->atoms[i].coor[j] - first->atoms[i].coor[j]);
      last->atoms[i].coor[j] = first->atoms[i].coor[j] + diff[i][j];
    }
  }

  sprintf(img, "img_%d.lammpstrj", 0);
  f0 = fopen(img, "w");
  LMPWriteFrame(f0, last, 0);
  fclose(f0);

  for (int i = 1; i < nimages; ++i)
  {
    sprintf(img, "img_%d.lammpstrj", i);
    f0 = fopen(img, "w");
    for (int j = 0; j < first->nat; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        last->atoms[j].coor[k] += diff[j][k];
      }
    }
    LMPWriteFrame(f0, last, 0);
    fclose(f0);
  }

  CrystalDelete(first);
  CrystalDelete(last);
  free(diff);

  return 0;
}
