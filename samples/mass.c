#include "fragments.h"
#include "lammpstrj.h"


int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat, pos, z; // , N, C, O, H;
  int spec[300];
  double dm[3], mass;
  for (int i = 0; i < 300; ++i)
  {
    printf("%d ", i );
  }

  // read the number of atoms and cell dimensions from lammpstrj
  LMPReadHeader(f, &t, &nat, dm);
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container
  Fragments *frg = FragmentsInit(bnd);

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    // read a frame
    LMPReadFrame(f, c);
    // assign the atoms to boxes
    BoxFill(c, box);
    // populate the bonding list
    BondingPopulate(c, box, bnd, checkBonding);
    // populate the fragments list
    FragmentsPopulate(bnd, frg);
    // Do the printing here!
    /* printf("%-8d\n", t); */

    pos = 0;
    // Iterate through all the fragments
    memset(spec, 0, 300*sizeof(int));
    for (int i = 0; i < frg->nfrags; ++i)
    {
      // set mass of fragment to zero
      mass = 0; // N = 0; C = 0; O = 0; H = 0;
      for (int j = 0; j < frg->lfrags[i]; ++j)
      {
        z = c->atoms[frg->frags[pos]].Z;
        switch (z) {
          case 1:
            /* H += 1; */
            mass +=  1.0079;
            break;
          case 6:
            /* C += 1; */
            mass += 12.0107;
            break;
          case 7:
            /* N += 1; */
            mass += 14.0067;
            break;
          case 8:
            /* O += 1; */
            mass += 15.9994;
            break;
        }
        pos += 1;
      }
      spec[(int)rint(mass)] += 1;
      /* printf("%8.3f ", mass); */
      /* if(C>0) printf("C"); */
      /* if(C>1) printf("%d", C); */
      /* if(H>0) printf("H"); */
      /* if(H>1) printf("%d", H); */
      /* if(O>0) printf("O"); */
      /* if(O>1) printf("%d", O); */
      /* if(N>0) printf("N"); */
      /* if(N>1) printf("%d", N); */
      /* printf("\n"); */
    }
    /* printf("\n\n"); */
    for (int i = 0; i < 300; ++i)
    {
      printf("%d ", spec[i]);
    }
    printf("\n");
    fflush(stdout);

    // read the header information for the next frame
    LMPReadHeader(f, &t, &nat, dm);
    // update cell size
    CrystalSetCell(c, dm);
    // Before moving to the next frame clear the bonding information
    BondingClear(bnd);
    // and the fragmentation information
    FragmentsClear(frg);
    // update the coarse graining box, in case cell dimensions changed
    BoxUpdate(c, box);
  }

  // free the memory
  fclose(f);
  FragmentsDelete(frg);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}

// vi: fdm=syntax
