#include "pdb.h"
#include "crystal.h"
#include <math.h>

// 10 / 4pi or 1/(4 pi dr)
#define NGRP  1             // number of atoms in the second group
#define MINR  2.0           // min dist to check gofr
#define NBINS 160           // number of bins
#define DR    0.1           // bin width
#define MAXR  NBINS*DR+MINR // DR*NBINS+MINR is the MAXR
#define RPI   0.0795774/DR  // 1 / (4 pi dr)
#define MOLS  300           // number of molecules
#define NAT   13            // number of atoms in a molecule

int main(int argc, char *argv[])
{

  FILE *pdb = fopen(argv[1], "r");

  int nat, nfr=0,
      bin,
      atm1,
      atm2,
      grp1   =   11  ,   // O
      grp2[] = { 11 }; // H2

  double dm[3],
         r,
         vol,
         gofr[NBINS],
         norm[NBINS];

  memset(gofr, 0, NBINS*sizeof(double));

  for (int i = 0; i < NBINS; ++i)
  {
    r = 1/(MINR+DR*(i+0.5)); // 1/radius
    norm[i] = RPI*r*r; // 1/4pir^2dr
  }

  PDBReadHeader(pdb, &nat, dm, 1);

  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);

  while(!feof(pdb))
  {
    fprintf(stderr, "Frame: %d\r", nfr++);
    vol = dm[0]*dm[1]*dm[2];
    PDBReadFrame(pdb, c);
    for (int i = 0; i < MOLS; ++i)
    {
      atm1 = i * NAT + grp1;
      for (int j = 0; j < MOLS; ++j)
      {
        if ( i == j ) continue;
        for (int k = 0; k < NGRP; ++k)
        {
          atm2 = j * NAT + grp2[k];
          r = CrystalDist(c, c->atoms[atm1].coor, c->atoms[atm2].coor);
          if ( r < MAXR )
          {
            bin = floor((r - MINR) / DR);
            gofr[bin] += vol * norm[bin];
          }
        }
      }
    }
    PDBReadHeader(pdb, &nat, dm, 0);
    CrystalSetCell(c, dm);
  }

  // free the memory
  fclose(pdb);
  CrystalDelete(c);

  for (int i = 0; i < NBINS; ++i)
  {
    printf("%5.2f %14.5f\n", MINR + DR * (i+0.5), gofr[i]/nfr/MOLS/MOLS/NGRP);
  }

  return 0;
}
