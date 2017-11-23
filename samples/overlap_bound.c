#include "bonding.h"
#include "xyz.h"
#include "lammpstrj.h"
#define MINDIST 0.9
#define MIN(X,Y) X<Y?X:Y;
#define MAX(X,Y) X>Y?X:Y;

/* check if there are any overlapping atoms in an xyz/lammpstrj file */

int CheckBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.2||(typ>2 && dist<1.4)||(typ>10 && dist<1.84))?1:0;
}

void PrintOut(Crystal *c, BondingInfo *bnd, int atm1, int atm2)
{
  double dist = CrystalDist(c, c->atoms[atm1].coor, c->atoms[atm2].coor);

  if(dist < MINDIST)
  {
    fprintf(stderr, "%4.2f\n", dist);

    fprintf(stderr, "%d %d ", c->atoms[atm1].Z, atm1);
    for (int i = 0; i < bnd->nbonds[atm1]; ++i)
    {
      fprintf(stderr, "%d ", bnd->bonds[atm1][i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "%d %d ", c->atoms[atm2].Z, atm2);
    for (int i = 0; i < bnd->nbonds[atm2]; ++i)
    {
      fprintf(stderr, "%d ", bnd->bonds[atm2][i]);
    }
    fprintf(stderr, "\n");
  }
}

int CheckOverlap(Crystal *c, CoarseBox *box)
{
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
  BondingInfo *bnd  = BondingInit(c);

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
    BondingPopulate(c, box, bnd, CheckBonding, 1);

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
            PrintOut(c, bnd, atm1, atm2);
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
              PrintOut(c, bnd, atm1, atm2);
            }
          }
        }
      }
    }

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
