// vi: fdm=syntax
#include "read.h"
#include "surface.h"

double beta(double coor[3])
{
  return sqrt(pow(coor[0]-50, 2) + pow(coor[1]-50, 2) + pow(coor[2]-50, 2));
}

int argmin(double arr[], int len)
{
  int j = 0;

  for (int i = 1; i < len; ++i)
  {
    if (arr[i] < arr[j])
    {
      j = i;
    }
  }
  return j;
}

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    printf("closest_surface traj.(xyz|lammpstrj) surfaces.in\n");
    return 1;
  }

  // read the surface.in
  // Format:
  // ax ay az // unit cell vectors
  // bx by bz
  // cx cy cz
  // nsurf  // number of surfaces
  // px1 py1 pz1 // point in the surface
  // vx1 vy1 vz1 // surface miller indices
  FILE *f = fopen(argv[2], "r");
  char *line = NULL;
  size_t len = 0;

  double cell[3][3];

  for (int i = 0; i < 3; ++i)
  {
    getline(&line, &len, f);
    sscanf(line, "%lf %lf %lf", &(cell[i][0]), &(cell[i][1]), &(cell[i][2]));
  }

  int nsurf;
  getline(&line, &len, f);
  sscanf(line, "%d", &nsurf);

  Surface s[nsurf];
  int miller[nsurf][3];

  for (int i = 0; i < nsurf; ++i)
  {
    getline(&line, &len, f);
    sscanf(line, "%lf %lf %lf", &(s[i].p[0]), &(s[i].p[1]), &(s[i].p[2]));
    getline(&line, &len, f);
    sscanf(line, "%d %d %d", &(miller[i][0]), &(miller[i][1]), &(miller[i][2]));
    miller_to_normal(miller[i], cell, s[i].n);
  }

  free(line);
  fclose(f);

  for (int i = 0; i < nsurf; ++i)
  {
    printf("{%7.2f %7.2f %7.2f} {%7.2f %7.2f %7.2f}\n",
        s[i].p[0], s[i].p[1], s[i].p[2],
        s[i].n[0], s[i].n[1], s[i].n[2]
        );
  }

  int t, xyz, nat, eof, counter;
  f = open_file(argv[1], &xyz);
  Crystal *c = read_file(f, xyz, &t);

  double dist[nsurf];
  int min, codes[10] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2};

  FILE *of = fopen("out.xyz", "w");
  do
  {
    fprintf(of, "%d\ncelldm 200 200 200\n", c->nat);
    for (int i = 0; i < c->nat; ++i)
    {
      for (int j = 0; j < nsurf; ++j)
      {
        dist[j] = fabs(SurfaceDist(s+j, c->atoms[i].coor));
      }
      min = argmin(dist, nsurf);
      fprintf(of, "%-2s %10.6f %10.6f %10.6f %10.6f %1d %10.6f\n",
          PT_Symbol(c->atoms[i].Z), c->atoms[i].coor[0], c->atoms[i].coor[1],
          c->atoms[i].coor[2], beta(c->atoms[i].coor), dist[min]<10?codes[min]:0, dist[min]);
      c->atoms[i].id = min+1;
    }

    for (int i = 0; i < 1; ++i)
    {
      eof = fgetc(f);
      if (eof == EOF) break;
      fseek(f, -1, SEEK_CUR);
      counter++;
      read_next(f, c, xyz, &t);
    }
  }
  while(eof != EOF);

  fclose(f);
  fclose(of);
  CrystalDelete(c);

  return 0;
}
