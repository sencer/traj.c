// vi: fdm=syntax
#include "lammpstrj.h"
#include "bonding.h"
#include "read.h"

int check_bonding(double len, int t1, int t2)
{
  int t = t1 + t2;
  return (
      len < 1.2 ||             // h-h
      (len < 1.4 && t > 2) ||  // o-h
      (len < 1.7 && t > 10) || // o-o, ti-h
      (len < 2.3 && t > 24) //|| // ti-o
      /* (len < 3.1 && t > 40)    // ti-ti */
      );
}

int main(int argc, char *argv[])
{
  int t, xyz;
  Crystal *c = read_single_frame(argv[1], &xyz, &t);
  CoarseBox *box = Box(c, 2.3);
  BondingInfo *bnd = Bonds(c, box, check_bonding, 1);
  BoxDelete(box);

  printf("set all [atomselect top all]\n");
  printf("$all set beta \"");
  for (int i = 0; i < c->nat; ++i)
  {
    printf("%7.4f ", sqrt( pow(c->atoms[i].coor[0]-50, 2) +
          pow(c->atoms[i].coor[1]-50, 2) + pow(c->atoms[i].coor[2]-50, 2) ));
  }
  printf("\"\n");
  CrystalDelete(c);

  printf("topo setbondlist \"");
  for (int i = 0; i < bnd->nat; ++i)
  {
    for (int j = 0; j < bnd->nbonds[i]; ++j)
    {
      if( i < bnd->bonds[i][j])
        printf(" {%d %d}", i, bnd->bonds[i][j]);
    }
  } 
  printf("\"\n");
  BondingDelete(bnd);
}
