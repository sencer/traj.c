#include "xyz.h"
#include "bonding.h"

#define MIN(X,Y) X<Y?X:Y;
#define MAX(X,Y) X>Y?X:Y;

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

  // read the number of atoms and cell dimensions from lammpstrj
  XYZReadHeader(f, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // read the file
  XYZReadFrame(f, c);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // assign the atoms to boxes
  BoxFill(c, box);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // populate the bonding list
  BondingPopulate(c, box, bnd, checkBonding, 1);
  /* int counter = 0;// bond_degree; */
  int atm0, atm1, bond_degree; 
  for (int i = 0; i < c->nat; ++i)
  {
    for (int j = 0; j < bnd->nbonds[i]; ++j)
    {
      atm0 = MIN(bnd->bonds[i][j], i);
      atm1 = MAX(bnd->bonds[i][j], i);
      if (bnd->nbonds[i] == 3 && bnd->nbonds[bnd->bonds[i][j]] == 3)
      {
        bond_degree = 2;
      }
      else
      {
        bond_degree = 1;
      }
      printf("%d %d %d\n", atm0+1, atm1+1, bond_degree);
    }
    /* if(c->atoms[i].Z == 6) */
    /* { */
    /*   if (bnd->nbonds[i] == 3) */
    /*   { */
    /*     c->atoms[i].id = 2; */
    /*   } */
    /*   else */
    /*   { */
    /*     for (int j = 0; j < 4; ++j) */
    /*     { */
    /*       if(c->atoms[bnd->bonds[i][j]].Z == 8) */
    /*       { */
    /*         c->atoms[i].id = c->atoms[bnd->bonds[i][j]].id; */
    /*       } */
    /*     } */
    /*   } */
    /* } */
  }
  /* FILE *out = fopen("out.xyz", "w"); */
  /* XYZWriteFrame(out, c); */
  /* fclose(out); */

  // free the memory
  fclose(f);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
