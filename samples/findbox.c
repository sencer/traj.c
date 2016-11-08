#include "lammpstrj.h"
#include "boxes.h"

int main(int argc, char *argv[])
{

  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat, comp[3];
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 4.205);

  // read a frame
  LMPReadFrame(f, c);
  // assign the atoms to boxes
  BoxFill(c, box);
  // Do the printing here!
  for (int i = 0; i < box->ntot; ++i)
  {
    if (box->binsn[i] == 0)
    {
      BoxGetComponents(box, i, comp);
      printf("O %.3f %.3f %.3f\n", (0.5+comp[0])*box->w[0]+0.32,
                                   (0.5+comp[1])*box->w[1]+0.32,
                                   (0.5+comp[2])*box->w[2]+0.32);
      printf("O %.3f %.3f %.3f\n", (0.5+comp[0])*box->w[0]-0.32,
                                   (0.5+comp[1])*box->w[1]-0.32,
                                   (0.5+comp[2])*box->w[2]-0.32);
    }
  }

  // free the memory
  fclose(f);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
