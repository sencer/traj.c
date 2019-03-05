// vi: fdm=syntax
#include <unistd.h>
#include "read.h"
#include "boxes.h"
#include "surface.h"
#include "hungarian.h"

#define CENTER 50.0
#define MAX_BETA 32.0
#define DO_TRANSFORM 0

void Transform(Crystal *c)
{
  double tmp[3], *coor;
  for (int i = 0; i < c->nat; ++i)
  {
    coor = c->atoms[i].coor;
    tmp[0] =  0.839*coor[0] + 0.405*coor[1] + 0.363*coor[2] - 45.549;
    tmp[1] =  0.000*coor[0] + 0.667*coor[1] - 0.745*coor[2] + 80.877;
    tmp[2] = -0.544*coor[0] + 0.625*coor[1] + 0.559*coor[2] + 26.955;
    memcpy(coor, tmp, 3*sizeof(double));
  }
}

double beta(double coor[3])
{
  double center[3] = {CENTER, CENTER, CENTER};

  return sqrt(pow(coor[0]-center[0], 2) +
      pow(coor[1]-center[1], 2) +
      pow(coor[2]-center[2], 2));
}

int argmin(double arr[], int len)
{
  int j = 0;

  for (int i = 1; i < len; ++i)
  {
    if (fabs(arr[i]) < fabs(arr[j]))
    {
      j = i;
    }
  }
  return j;
}

int MapAtoms(Crystal *ref, CoarseBox *rbox, Crystal *tar, CoarseBox *tbox)
{
  int (*costs)[ref->nat] = malloc(ref->nat * sizeof(*costs)),
      boxlist[ref->nat], irbox, itbox, ir, it, neigh[27], nneigh, mapped=0;

  BoxesOfAtoms(ref, rbox, boxlist);

  // fill the cost matrix with a high cost initially
  for (int i = 0; i < ref->nat; ++i)
  {
    for (int j = 0; j < ref->nat; ++j)
    {
      costs[i][j] = 99999;
    }
  }

  // build the real cost matrix
  for (ir = 0; ir < ref->nat; ++ir)
  {
    irbox = boxlist[ir];
    nneigh = GetNeighboringBoxes(rbox, irbox, neigh, 0);

    for (int j = 0; j < nneigh; ++j)
    {
      itbox = neigh[j];

      for (int k = 0; k < tbox->binsn[itbox]; ++k)
      {
        it = tbox->bins[itbox][k];

        if (ref->atoms[ir].Z == tar->atoms[it].Z &&
            beta(tar->atoms[it].coor) < MAX_BETA)
        {
          costs[ir][it] = (int)(round(100*sqrt(
                  pow(ref->atoms[ir].coor[0] - tar->atoms[it].coor[0], 2) +
                  pow(ref->atoms[ir].coor[1] - tar->atoms[it].coor[1], 2) +
                  pow(ref->atoms[ir].coor[2] - tar->atoms[it].coor[2], 2)
                  )));
        }
      }
    }
  }

  // solve the bipartite matching problem
  hungarian_problem_t p;
  int matrix_size = hungarian_init(&p, ref->nat, ref->nat, costs, 0);
  hungarian_solve(&p);

  for (ir = 0; ir < ref->nat; ++ir)
  {
    for (it = 0; it < tar->nat; ++it)
    {
      if ( p.assignment[ir][it] == 1)
      {
        ref->atoms[ir].id = it;
        tar->atoms[it].id = costs[ir][it];
      }
    }
  }

  hungarian_free(&p);
  free(costs);

  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    printf("np_rmsd.c traj.xyz surfaces.in\n");
    return 1;
  }

  if (access(argv[1], F_OK) == -1)
  {
    printf("%s does not exist.\n", argv[1]);
    return 1;
  }

  if (access(argv[2], F_OK) == -1)
  {
    printf("%s does not exist.\n", argv[2]);
    return 1;
  }

  // read the surface definitions from argv[2]
  FILE *fsurf = fopen(argv[2], "r");
  char *line = NULL;
  size_t len = 0;

  // read unit cell vectors
  double cell[3][3];

  for (int i = 0; i < 3; ++i)
  {
    getline(&line, &len, fsurf);
    sscanf(line, "%lf %lf %lf", cell[i], cell[i]+1, cell[i]+2);
  }

  // read number of surface definitions
  int nsurf;
  getline(&line, &len, fsurf);
  sscanf(line, "%d", &nsurf);

  // read the surface definitions
  int miller[nsurf][3], equiv[nsurf];
  Surface s[nsurf];

  for (int i = 0; i < nsurf; ++i)
  {
    getline(&line, &len, fsurf);
    sscanf(line, "%d %lf %lf %lf %d %d %d", equiv+i, s[i].p, s[i].p+1,
        s[i].p+2, miller[i], miller[i]+1, miller[i]+2);
    // use the unit cell definition to convert miller indices
    // to surface normal in real (Angstrom) units
    miller_to_normal(miller[i], cell, s[i].n);
  }

  free(line);
  fclose(fsurf);

  // start reading trajectory
  int t, xyz, nat, eof, counter = 1;

  FILE *ftraj = open_file(argv[1], &xyz),
       *fout  = fopen("out1.xyz", "w");

  Crystal *ref = read_file(ftraj, xyz, &t);
#if DO_TRANSFORM
  Transform(ref);
#endif
  CoarseBox *rbox = Box(ref, 3.0);

  // some data holders for choosing the closest surface for each atom
  double sdist[nsurf];
  int isurf1, isurf2, flag;

  // create data holders for target frames
  Crystal *tar = CrystalInit(ref->nat, ref->dm);
  CoarseBox *tbox = BoxInit(tar, 3.0);

  do
  {
    for (int i = 0; i < 1; ++i)
    {
      eof = fgetc(ftraj);
      if (eof == EOF) break;
      fseek(ftraj, -1, SEEK_CUR);
      read_next(ftraj, tar, xyz, &t);
#if DO_TRANSFORM
      Transform(tar);
#endif
      counter++;
      printf("\rProcessing frame %d", counter);
      fflush(stdout);
    }
    BoxClear(tbox);
    BoxFill(tar, tbox);

    MapAtoms(ref, rbox, tar, tbox);

    fprintf(fout, "0\n");

    Atom at;
    for (int ir = 0; ir < tar->nat; ++ir)
    {
      at = tar->atoms[ref->atoms[ir].id];
      flag = 0;

      for (int j = 0; j < nsurf; ++j)
      {
        sdist[j] = SurfaceDist(s+j, at.coor);
        if (sdist[j] > 1)
        {
          flag = 1;
          break;
        }
      }

      if (flag) continue;

      isurf1 = argmin(sdist, nsurf-2);   // {101} facets
      isurf2 = nsurf - 2 + argmin(sdist+nsurf-2, 2); // {001} facets

      fprintf(fout, "%10.6f %10.6f %10.6f\n", sdist[isurf1], sdist[isurf2],
          at.id/100.0);
    }
  }
  while (eof!=EOF);

  printf("\n");

  CrystalDelete(ref);
  CrystalDelete(tar);
  BoxDelete(rbox);
  BoxDelete(tbox);
  fclose(ftraj);
  fclose(fout);

  return 0;
}
