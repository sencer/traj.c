#include "lammpstrj.h"
#include "fragments.h"
#define SHEET 100

/*
 * Reads a trajectory to analyze the functional groups on FGSs, and the small
 * molecules formed from these functional groups.
 */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0 ||
         (typ>2 && dist<1.4) ||
         (typ>10 && dist<1.9))?1:0; // ||
         /* (t1==6 && t2==6 && dist<2.1))?1:0; */
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat;
  double dm[3];

  LMPReadHeader(f, &t, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  CoarseBox *box = BoxInit(c, 2.0);
  BondingInfo *bnd = BondingInit(c);
  FragmentsInfo *frg = FragmentsInit(bnd);

  printf("# Time    C/O    H/O     fOH  fEpox fEthr fCarb fOthr  OH    H2O   H2O2  O2    CO    CO2  Cxx Other TotMols Vtot\n");

  int pos, id, isBonded, n_cc_bonds, n_removed_c, n_c_sheet,
      elem[3],       // H, O, C
      ratio[3],      // H, O, C
      functional[5], // OH, Epoxy, In-plane ether, Carbonyl, Other
      molecules[8];  // OH, H2O, H2O2, O2, CO, CO2, other_c, other
  Atom *atm;
  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    // read a frame
    LMPReadFrame(f, c);
    BoxFill(c, box);
    BondingPopulate(c, box, bnd, checkBonding, 1);
    FragmentsPopulate(bnd, frg);

    memset(ratio, 0, 3*sizeof(int));
    memset(functional, 0, 5*sizeof(int));
    memset(molecules, 0, 8*sizeof(int));
    pos = 0;
    n_c_sheet = 0;
    n_cc_bonds = 0;
    n_removed_c = 0;

    /* for (int i = 0; i < c->nat; ++i) */
    /* { */
    /*   if (c->atoms[i].Z == 6) */
    /*   { */
    /*     for (int j = 0; j < bnd->nbonds[i]; ++j) */
    /*     { */
    /*       if (c->atoms[bnd->bonds[i][j]].Z == 6) */
    /*       { */
    /*         n_cc_bonds++; */
    /*       } */
    /*     } */
    /*   } */
    /* } */

    for (int i = 0; i < frg->nfrags; ++i)
    {
      memset(elem, 0, 3*sizeof(int));
      for (int j = 0; j < frg->lfrags[i]; ++j)
      {
        id = frg->frags[pos++];
        /* if ( id >= 471) continue; */
        atm = &(c->atoms[id]);
        switch (atm->Z)
        {
          case 1 : elem[0]++; break;
          case 8 : elem[1]++; break;
          case 6 : elem[2]++; break;
          default: printf("Unexpected atom.\n");
        }
        // Check the functional groups.
        if(frg->lfrags[i]>SHEET && atm->Z == 8)
        {
          switch (bnd->nbonds[id])
          {
            case 1:
              if(c->atoms[bnd->bonds[id][0]].Z==6)
              {
                functional[3]++;
              }
              else
              {
                functional[4]++;
              }
              break;
            case 2:
              if( c->atoms[bnd->bonds[id][0]].Z
                + c->atoms[bnd->bonds[id][1]].Z == 12 )
              {
                isBonded = 0;
                for (int k = 0; k < bnd->nbonds[bnd->bonds[id][0]]; ++k)
                {
                  if(bnd->bonds[bnd->bonds[id][0]][k] == bnd->bonds[id][1])
                  {
                    isBonded = 1;
                  }
                }
                if (isBonded)
                {
                  functional[1]++;
                }
                else
                {
                  functional[2]++;
                }
              }
              else if( c->atoms[bnd->bonds[id][0]].Z
                     + c->atoms[bnd->bonds[id][1]].Z == 7 )
              {
                functional[0]++;
              }
              else
              {
                functional[4]++;
              }
              break;
            default:
              functional[4]++;
          }
        }
      }
      if(frg->lfrags[i] > SHEET)
      {
        ratio[0] += elem[0];
        ratio[1] += elem[1];
        ratio[2] += elem[2];

        pos = 0;
        for (int j = 0; j < frg->lfrags[i]; ++j)
        {
          id = frg->frags[pos++];
          /* if ( id >= 471) continue; */
          atm = &(c->atoms[id]);
          if (atm->Z == 6)
          {
            n_c_sheet++;
            for (int k = 0; k < bnd->nbonds[id]; ++k)
            {
              if (bnd->bonds[id][k]<id && c->atoms[bnd->bonds[id][k]].Z == 6)
              {
                n_cc_bonds++;
              }
            }
          }
        }

      }
      else
      {
             if( elem[0]==1 && elem[1]==1 && elem[2]==0 ) molecules[0]++;
        else if( elem[0]==2 && elem[1]==1 && elem[2]==0 ) molecules[1]++;
        else if( elem[0]==2 && elem[1]==2 && elem[2]==0 ) molecules[2]++;
        else if( elem[0]==0 && elem[1]==2 && elem[2]==0 ) molecules[3]++;
        else if( elem[0]==0 && elem[1]==1 && elem[2]==1 ) molecules[4]++;
        else if( elem[0]==0 && elem[1]==2 && elem[2]==1 ) molecules[5]++;
        else if(                             elem[2]> 0 ) molecules[6]++;
        else                                              molecules[7]++;
        n_removed_c += elem[2];
      }
    }

    printf("%-8d %5.3f %5.3f %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %6d %10.3f %4d %4d %5d\n",
        t, 1.0*ratio[2]/ratio[1], 1.0*ratio[0]/ratio[1], functional[0],
        functional[1], functional[2], functional[3], functional[4],
        molecules[0], molecules[1], molecules[2], molecules[3], molecules[4],
        molecules[5], molecules[6] ,molecules[7], frg->nfrags,
        c->dm[0]*c->dm[1]*c->dm[2], n_cc_bonds, n_removed_c, n_c_sheet);

    fflush(stdout);

    LMPReadHeader(f, &t, &nat, dm);
    if (nat != c->nat)
    {
      CrystalDelete(c);
      BondingDelete(bnd);
      FragmentsDelete(frg);
      c = CrystalInit(nat, dm);
      bnd = BondingInit(c);
      frg = FragmentsInit(bnd);
    }
    else
    {
      CrystalSetCell(c, dm);
      BondingClear(bnd);
      FragmentsClear(frg);
    }
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
