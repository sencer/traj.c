#define _GNU_SOURCE
#include "lammpstrj.h"

int LMPReadHeader(FILE *f, int *t, int *nat, double *dm)
{
  char *line = NULL;
  double tmp1, tmp2;
  size_t len = 0;
  // TODO Check the section headers being skipped now

  getline(&line, &len, f);
  getline(&line, &len, f);
  sscanf(line, "%d", t);

  getline(&line, &len, f);
  getline(&line, &len, f);
  sscanf(line, "%d", nat);

  getline(&line, &len, f);
  for (int i = 0; i < 3; ++i)
  {
    getline(&line, &len, f);
    sscanf(line, "%lf %lf", &tmp1, &tmp2);
    dm[i] = tmp2 - tmp1;
  }

  free(line);
  return 0;
}

int LMPReadFrame(FILE *f, Crystal *c)
{
  char *line = NULL;
  char str[5];
  unsigned char sym[3];
  size_t len = 0;

  getline(&line, &len, f);
  for (int i = 0; i < c->nat; ++i)
  {
    getline(&line, &len, f);
    sscanf(line, "%s %lf %lf %lf", str, c->atoms[i].coor,
                                        c->atoms[i].coor + 1,
                                        c->atoms[i].coor + 2);

    sscanf(str, "%[^0-9]%d", sym, &(c->atoms[i].id));

    for (int j = 0; j < 3; ++j)
    {
      c->atoms[i].coor[j] -= floor(c->atoms[i].coor[j]/c->dm[j]) * c->dm[j];
    }
    c->atoms[i].Z = PT_AtomicNumber(sym);
  }

  free(line);
  return 0;
}

int LMPWriteFrame(FILE *f, Crystal *c, int t)
{
  fprintf(f, "ITEM: TIMESTEP\n%d\n", t);
  fprintf(f, "ITEM: NUMBER OF ATOMS\n%d\n", c->nat);
  fprintf(f, "ITEM: BOX BOUNDS pp pp pp\n");
  for (int i = 0; i < 3; ++i)
  {
    fprintf(f, "0.0 %12.8f\n", c->dm[i]);
  }
  fprintf(f, "ITEM: ATOMS element x y z\n");
  for (int i = 0; i < c->nat; ++i)
  {

    if (c->atoms[i].id > 0)
    {
      fprintf(f, "%2s%-3d", PT_Symbol(c->atoms[i].Z), c->atoms[i].id);
    }
    else
    {
      fprintf(f, " %-4s", PT_Symbol(c->atoms[i].Z));
    }

    fprintf(f, " %12.9f %12.9f %12.9f\n", c->atoms[i].coor[0],
                                          c->atoms[i].coor[1],
                                          c->atoms[i].coor[2]);
  }

  return 0;
}

// TODO two alternatives, one wrapped, one unwrapped coordinates
