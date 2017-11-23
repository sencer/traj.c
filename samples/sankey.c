#include "lammpstrj.h"
#include "fragments.h"
#define SHEET 50

/*
 * Reads a frame and assigns one of these to O atoms in the order they are 
 * encountered:
 * 0 -> Epoxy
 * 1 -> OH
 * 2 -> Ether
 * 3 -> Carbonyl
 * 4 -> Other surf
 * 5 -> H2O
 * 6 -> O2 / H2O2
 * 7 -> CO / CO2
 * 8 -> other mols
 */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0 ||
         (typ>2 && dist<1.4) ||
         (typ>10 && dist<1.9)?1:0||
         (t1==6 && t2==6 && dist<2.1))?1:0;
}

int isBound(BondingInfo *bnd, int a1, int a2)
{
  int bound = 0;
  for (int i = 0; i < bnd->nbonds[a1]; ++i)
  {
    if (bnd->bonds[a1][i] == a2)
    {
      bound = 1;
      break;
    }
  }
  return bound;
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat;
  double dm[3];

  LMPReadHeader(f, &t, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container
  FragmentsInfo *frg = FragmentsInit(bnd);

  LMPReadFrame(f, c);
  BoxFill(c, box);
  BondingPopulate(c, box, bnd, checkBonding, 1);
  FragmentsPopulate(bnd, frg);

  int pos = 0;
  // first assign each atom to its fragment
  for (int i = 0; i < frg->nfrags; ++i)
  {
    for (int j = 0; j < frg->lfrags[i]; ++j)
    {
      c->atoms[frg->frags[pos++]].id = i;
    }
  }

  int tmp, elem[3], id;
  for (int i = 0; i < c->nat; ++i)
  {
    // only look at O atoms
    if (c->atoms[i].Z != 8) continue;
    id = c->atoms[i].id;

    // if the atom is on the FGS
    if (frg->lfrags[id] > SHEET)
    {
      switch (bnd->nbonds[i]) {
        case 1:
          if(c->atoms[bnd->bonds[i][0]].Z == 6)
          {
            printf("3\n"); // carbonyl
          }
          else
          {
            /* fprintf(stderr, "Warning 4: atom %d, fragment %d, nbonds 1.\n", i, id); */
            printf("4\n"); // other surface species
          }
          break;
        case 2:
          tmp = c->atoms[bnd->bonds[i][0]].Z + c->atoms[bnd->bonds[i][1]].Z;
          if (tmp == 7)
          {
            printf("1\n"); // OH
          }
          else if (tmp == 12)
          {
            if(isBound(bnd, bnd->bonds[i][0], bnd->bonds[i][1]))
            {
              printf("0\n"); //epoxy
            }
            else
            {
              printf("2\n"); // ether
            }
          }
          else
          {
            /* fprintf(stderr, "Warning 4: atom %d, fragment %d, nbonds 2.\n", i, id); */
            printf("4\n"); // other surf
          }
          break;
        default:
          printf("4\n"); // other surf
      }
    }
    // if it is a small molecule
    else
    {
      memset(elem, 0, 3*sizeof(int));

      // find the beginning of the fragment
      pos = 0;
      for (int j = 0; j < id; ++j)
      {
        pos += frg->lfrags[j];
      }
      // walk the atoms in the fragment
      for (int j = 0; j < frg->lfrags[id]; ++j)
      {
        /*
         *  1%3 -> 1 == H
         *  6%3 -> 0 == C
         *  8%3 -> 2 == O
         */
        elem[c->atoms[frg->frags[pos+j]].Z%3]++;
      }

      if ( elem[0]==0 && (elem[1]==1||elem[1]==2) && (elem[2]==2||elem[2]==1) ) // H2O/OH+H2O2/HO2
      {
        printf("5\n");
      }
      else if ( elem[0]==0 && elem[1]==0 && (elem[2]==3||elem[2]==2) ) // O2/O3
      {
        printf("6\n");
      }
      else if ( elem[0]==1 && elem[1]==0 && (elem[2]==1||elem[2]==2) ) // CO, CO2
      {
        /* if(frg->lfrags[id] == 3) */
        /*   fprintf(stderr, "%d %d %d\n", frg->frags[pos], frg->frags[pos+1], frg->frags[pos+2]); */
        /* else */
        /* fprintf(stderr, "%d %d\n", frg->frags[pos], frg->frags[pos+1]); */
        printf("7\n");
      }
      else
      {
        printf("8\n"); // other mols
      }
    }
  }


  // free the memory
  fclose(f);
  FragmentsDelete(frg);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
