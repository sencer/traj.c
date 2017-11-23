#include "xyz.h"
#include "crystal.h"

/*
 * Reads a trajectory and calculates the diffusion of each atom.
 */

void SwapAtoms(Crystal *c1, Crystal *c2)
{
  void* tmp = c1->atoms;
  c1->atoms = c2->atoms;
  c2->atoms = tmp;
}

double dist(double a[3], double b[3])
{
  return sqrt( pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2));
}

int main(int argc, char *argv[])
{
  // trajectory to read
  FILE *f = fopen(argv[1], "r");

  int nat, // number of atoms
     *irrelevant; // an array to eliminate atoms that diffused out of limits

  double dm[3],       // cell dimensions
        *diff,        // total length of diffusion
         center[3] = {50,50,50};

  XYZReadHeader(f, &nat, dm);

  irrelevant  = malloc(nat * sizeof(int));
  memset(irrelevant, 0, nat*sizeof(int));

  diff  = malloc(nat * sizeof(double));
  memset(diff, 0, nat*sizeof(double));

  Crystal *cur = CrystalInit(nat, dm);
  Crystal *prev = CrystalInit(nat, dm);
  XYZReadFrame(f, prev);

  while(!feof(f))
  {
    // read the following frame
    XYZReadHeader(f, &nat, dm);
    XYZReadFrame(f, cur);

    // for each atom
    for (int i = 0; i < nat; ++i)
    {
      diff[i] += dist(cur->atoms[i].coor, prev->atoms[i].coor);
      if ( dist(cur->atoms[i].coor, center) > 33.0 )
      {
        irrelevant[i] = 1;
      }
    }
    // swap the pointers for current and previous frames
    SwapAtoms(prev, cur);
  }
  for (int i = 0; i < nat; ++i)
  {
    if(!irrelevant[i])
    {
      printf("%4d %6.3f\n", i, diff[i]);
    }
  }


  free(diff);
  fclose(f);
  CrystalDelete(prev);
  CrystalDelete(cur);

  return 0;
}
