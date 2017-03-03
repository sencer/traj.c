#include "bonding.h"
#include "xyz.h"

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t = 0, nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  XYZReadHeader(f, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container

  // now, while we didn't hit the end of file
  int flag;
  while(!feof(f))
  {
    // read a frame
    XYZReadFrame(f, c);
    // assign the atoms to boxes
    BoxFill(c, box);
    // populate the bonding list
    BondingPopulate(c, box, bnd, checkBonding);
    // populate the fragments list
    // Do the printing here!
    printf("%d ", t);
    flag = 0;
    for (int i = 144; i < 216; ++i)
    {
      if (bnd->nbonds[i] == 3 &&
          bnd->bonds[i][0]>215 &&
          bnd->bonds[i][1]>215 &&
          bnd->bonds[i][2]>215
          )
      {
        printf("%d %d %d %d ", i,
            bnd->bonds[i][0],
            bnd->bonds[i][1],
            bnd->bonds[i][2]
            );
        flag = 1;
      }
    }
    if(!flag) printf("10000 ");
    printf("\n");

    // read the header information for the next frame
    XYZReadHeader(f, &nat, dm);
    // update cell size
    CrystalSetCell(c, dm);
    // Before moving to the next frame clear the bonding information
    BondingClear(bnd);
    // and the fragmentation information
    // update the coarse graining box, in case cell dimensions changed
    BoxUpdate(c, box);
    t++;
  }

  // free the memory
  fclose(f);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
