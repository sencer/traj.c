#include "overlap.h"
#define MINDIST 0.65

int CheckOverlap(Crystal *c, CoarseBox *box)
{
  int icomp[3], jcomp[3], jbox, atm1, atm2;
  double dist;

  // for each box
  for (int ibox = 0; ibox < box->ntot; ++ibox)
  {
    // if there are atoms in the box
    if (box->binsn[ibox])
    {
      // cross-check bonding between the atoms in box ibox
      for (int i = 0; i < box->binsn[ibox]-1; ++i)
      {
        atm1 = box->bins[ibox][i];
        for (int j = i+1; j < box->binsn[ibox]; ++j)
        {
          atm2 = box->bins[ibox][j];
          dist = CrystalDist(c, atm1, atm2);
          if(dist < MINDIST)
          {
            printf("%d and %d overlapping. (%f)\n", atm1, atm2, dist);
          }
        }
      }
      // now check the neighboring boxes - only the 13 of them in (+) direction
      BoxGetComponents(box, ibox, icomp);
      for (int j = 14; j < 27; ++j)
      {
        jcomp[0] = icomp[0] + j / 9 - 1;
        jcomp[1] = icomp[1] + (j % 9) / 3 - 1;
        jcomp[2] = icomp[2] + j % 3 - 1;
        jbox = BoxGetIndice(box, jcomp);
        // for each atom in box jbox check bonding
        for (int k = 0; k < box->binsn[jbox]; ++k)
        {
          atm1 = box->bins[jbox][k];
          for (int i = 0; i < box->binsn[ibox]; ++i)
          {
            atm2 = box->bins[ibox][i];
            dist = CrystalDist(c, atm1, atm2);
            if(dist < MINDIST)
            {
              printf("%d and %d overlapping. (%f)\n", atm1, atm2, dist);
            }
          }
        }
      }
    }
  }
  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    printf("Please give a file name.\n");
    return -1;
  }
  FILE *f = fopen(argv[1], "r");
  int nat, t;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  /* XYZReadHeader(f, &nat, dm); */
  LMPReadHeader(f, &t, &nat, dm);
  printf("%f %f %f\n", dm[0], dm[1], dm[2]);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  CoarseBox *box = BoxInit(c);
  // read a frame
  /* XYZReadFrame(f, c); */
  while(!feof(f))
  {
    LMPReadFrame(f, c);
    // ...and a coarse graining box
    // assign the atoms to boxes
    BoxFill(c, box);
    CheckOverlap(c, box);
    LMPReadHeader(f, &t, &nat, dm);
    CrystalSetCell(c, dm);
    BoxUpdate(c, box);
  }

  // free the memory
  fclose(f);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
