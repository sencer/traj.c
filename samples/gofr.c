#define _GNU_SOURCE
#include "read.h"
#include <math.h>

/* SEE ALSO THE CONSTANTS BELOW
 * Calculates g(r) for different species in a trajectory. Requires a 'gofr.in' 
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
 */

// 10 / 4pi or 1/(4 pi dr)
#define MINR  0.5           // min dist to check gofr
#define NBINS 180           // number of bins
#define DR    0.1           // bin width
#define MAXR  NBINS*DR+MINR // DR*NBINS+MINR is the MAXR

double volume(double r)
{
  return 4.188666*pow(r, 3);
}

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
  int ngofrs;
  fscanf(in, "%d\n", &ngofrs);

  int eof, bin, xyz, t, nf, ngrp1, ngrp2, *grp1, *grp2;
  double norm[NBINS], gofr[NBINS], r, vol;
  FILE *f = open_file(argv[1], &xyz), *of;
  Crystal *c = read_file(f, xyz, &t);

  norm[0] = volume(MINR + DR*0.5);
  for (int i = 1; i < NBINS; ++i)
  {
    norm[i] = volume(MINR + (i+0.5)*DR) - norm[i-1];
  }

  for (int i = 0; i < NBINS; ++i)
  {
    norm[i] = 1 / norm[i];
  }


  for (int icalc = 0; icalc < ngofrs; ++icalc)
  {
    fscanf(in, "%s\n", outfile);
    fscanf(in, "%d %d", &ngrp1, &ngrp2);
    grp1 = malloc(ngrp1 * sizeof(int));
    read_group(in, ngrp1, grp1);
    grp2 = malloc(ngrp2 * sizeof(int));
    read_group(in, ngrp2, grp2);

    memset(gofr, 0, NBINS*sizeof(double));
    fseek(f, 0, SEEK_SET);
    eof = 0;
    nf = 0;

    while(eof != EOF)
    {
      read_next(f, c, xyz, &t);
      fprintf(stderr, "File %s, Frame %d\r", outfile, nf++);
      vol = c->dm[0] * c->dm[1] * c->dm[2];
      for (int i = 0; i < ngrp1; ++i)
      {
        for (int j = 0; j < ngrp2; ++j)
        {
          r = CrystalDist(c, c->atoms[grp1[i]].coor, c->atoms[grp2[j]].coor);
          if (r < MAXR)
          {
            bin = floor((r - MINR) / DR);
            gofr[bin] += vol * norm[bin];
          }
        }
      }
      eof = fgetc(f);
      fseek(f, -1, SEEK_CUR);
    }
    of = fopen(outfile, "w");
    for (int i = 0; i < NBINS; ++i)
    {
      fprintf(of, "%5.2f %14.5f\n", MINR + DR*(i+0.5),
                                    gofr[i]/nf/ngrp1/ngrp2);
    }
    fclose(of);
    free(grp1);
    free(grp2);
  }

  fclose(f);
  fclose(in);
  CrystalDelete(c);

  return 0;
}
