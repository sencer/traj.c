#define _GNU_SOURCE
#include "pdb.h"
#include "crystal.h"
#include <math.h>

/* Calculates g(r) for different species in a trajectory. Requires a 'gofr.in' 
 * to decide on the species.
 * Format of gofr.in is as follows:
 * 2                -> number of g(r) calculations to do.
 * filename.1       -> file name to save the results of the first calculation
 * 1 1              -> number of atoms in the first and second group
 * 11               -> index of the atoms in the first group
 * 11               -> index of the atoms in the seconds group
 * filename.2       -> same for the second group.
 * ...
 *
 * Currently the implementation is terrible - but at least it works.
 * Reads the trajectory again and again for each calculation.
 * But can be improved easily.
 * SEE ALSO THE CONSTANTS BELOW
 */

// 10 / 4pi or 1/(4 pi dr)
#define MINR  0.0           // min dist to check gofr
#define NBINS 180           // number of bins
#define DR    0.1           // bin width
#define MAXR  NBINS*DR+MINR // DR*NBINS+MINR is the MAXR
#define RPI   0.0795774/DR  // 1 / (4 pi dr)
#define MOLS  300           // number of molecules
#define NAT   13            // number of atoms in a molecule

int read_group(FILE *f, int len, int *arr)
{
  for (int i = 0; i < len; ++i)
  {
    fscanf(f, "%d", arr+i);
  }

  return 0;
}

int main(int argc, char *argv[])
{

  FILE *in = fopen("gofr.in", "r");

  char outfile[10];
  int n_gofr_calculations;
  fscanf(in, "%d\n", &n_gofr_calculations);


  FILE *pdb, *f;
  Crystal *c;

  int nat, nfr,
      bin,
      atm1,
      atm2,
      ngrp1,
      ngrp2,
      *grp1,
      *grp2;

  double dm[3], r, vol, gofr[NBINS], norm[NBINS];

  for (int i = 0; i < NBINS; ++i)
  {
    r = 1/(MINR+DR*(i+0.5)); // 1/radius
    norm[i] = RPI*r*r; // 1/4pir^2dr
  }


  // setup the calculation.
  for (int ncalc = 0; ncalc < n_gofr_calculations; ++ncalc)
  {
    fscanf(in, "%s\n", outfile);
    fscanf(in, "%d %d\n", &ngrp1, &ngrp2);
    grp1 = malloc(ngrp1 * sizeof(int));
    read_group(in, ngrp1, grp1);
    grp2 = malloc(ngrp2 * sizeof(int));
    read_group(in, ngrp2, grp2);


    // a stupid implementation which reads the trajectory again and again.
    pdb = fopen(argv[1], "r");
    PDBReadHeader(pdb, &nat, dm, 1);

    // initialize a crystal structure with these dimensions
    c = CrystalInit(nat, dm);

    nfr = 0;
    memset(gofr, 0, NBINS*sizeof(double));
    fprintf(stderr, "                                                      \r");
    while(!feof(pdb))
    {
      fprintf(stderr, "File %s, Frame: %d\r", outfile, nfr++);
      vol = dm[0]*dm[1]*dm[2];
      PDBReadFrame(pdb, c);
      for (int i = 0; i < MOLS; ++i)
      {
        for (int j = 0; j < MOLS; ++j)
        {
          if ( i == j ) continue;
          for (int l = 0; l < ngrp1; ++l)
          {
            atm1 = i * NAT + grp1[l];
            for (int k = 0; k < ngrp2; ++k)
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
      }
      PDBReadHeader(pdb, &nat, dm, 0);
      CrystalSetCell(c, dm);
    }

    f = fopen(outfile, "w");
    for (int i = 0; i < NBINS; ++i)
    {
      fprintf(f, "%5.2f %14.5f\n", MINR + DR * (i+0.5), gofr[i]/nfr/MOLS/MOLS/ngrp1/ngrp2);
    }
    // free the memory
    fclose(pdb);
    fclose(f);
    CrystalDelete(c);
    free(grp1);
    free(grp2);
  }

  fclose(in);
  return 0;
}
