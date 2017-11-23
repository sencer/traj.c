#include "boxes.h"
#include "xyz.h"
#include "lammpstrj.h"
#define MINDIST 0.5
#define MIN(X,Y) X<Y?X:Y;
#define MAX(X,Y) X>Y?X:Y;

/* check if there are any overlapping atoms in an xyz/lammpstrj file */

void PrintOut(Crystal *c, int atm1, int atm2)
{
  double dist = CrystalDist(c, c->atoms[atm1].coor, c->atoms[atm2].coor);
  int small = MIN(atm1, atm2);

  if ((small == 1 && dist<0.7)||dist<1.3)
  {
    /* printf("%5.2f set vsel [atomselect top \"same fragment as index %d %d\"]\n",dist , atm1, atm2); */
    printf("%5.2f set vsel [atomselect top \"index %d %d\"]\n", dist , atm1, atm2);
  }
  /* if(dist < MINDIST && c->atoms[atm1].Z == 1 && c->atoms[atm2].Z == 1) */
  /* { */
  /* /1*   /2* if(dist < 0.1) *2/ *1/ */
  /* /1*   /2* { *2/ *1/ */
  /* /1*   /2*   fprintf(stderr, "0.xyz:%d\n", atm2+3); *2/ *1/ */
  /* /1*   /2* } *2/ *1/ */
  /*   fprintf(stderr, "%d %d\n", atm1, atm2); */
  /* } */
}

int CheckOverlap(Crystal *c, CoarseBox *box)
{
  int icomp[3], jcomp[3], jbox, atm1, atm2;

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
          PrintOut(c, atm1, atm2);
        }
      }
      // now check the neighboring boxes - only the 13 of them in (+) direction
      BoxGetComponents(box, ibox, icomp);
      for (int j = 14; j < 27; ++j)
      {
        jcomp[0] = icomp[0] + j / 9 - 1;
        jcomp[1] = icomp[1] + (j % 9) / 3 - 1;
        jcomp[2] = icomp[2] + j % 3 - 1;
        jbox = BoxGetIndice(box, jcomp, 1);
        // for each atom in box jbox check bonding
        for (int k = 0; k < box->binsn[jbox]; ++k)
        {
          atm1 = box->bins[jbox][k];
          for (int i = 0; i < box->binsn[ibox]; ++i)
          {
            atm2 = box->bins[ibox][i];
            PrintOut(c, atm1, atm2);
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

  int xyz = strcmp(strrchr(argv[1], '.'), ".xyz") != -1;
  // read the number of atoms and cell dimensions from lammpstrj
  if (xyz)
    XYZReadHeader(f, &nat, dm);
  else
    LMPReadHeader(f, &t, &nat, dm);

  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  CoarseBox *box = BoxInit(c, 2.0);

  // will end with an error for lammpstrj files, since EOF in LMPReadHeader
  while(!feof(f))
  {
    // read a frame
    if (xyz)
      XYZReadFrame(f, c);
    else
      LMPReadFrame(f, c);
    // ...and a coarse graining box
    // assign the atoms to boxes
    BoxFill(c, box);
    CheckOverlap(c, box);

    if (xyz)
      XYZReadHeader(f, &nat, dm);
    else
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
