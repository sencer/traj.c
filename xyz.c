#define _GNU_SOURCE
#include "xyz.h"

int XYZReadHeader(FILE *f, int *nat, double *dm)
{
  char *line = NULL;
  size_t len = 0;
  // TODO Check the section headers being skipped now

  getline(&line, &len, f);
  sscanf(line, "%d", nat);
  getline(&line, &len, f);
  sscanf(line, "%*s %lf %lf %lf %*f %*f %*f", dm, dm+1, dm+2);

  free(line);
  return 0;
}

int XYZReadFrame(FILE *f, Crystal *c)
{
  char *line = NULL;
  char str[3];
  size_t len = 0;

  for (int i = 0; i < c->nat; ++i)
  {
    getline(&line, &len, f);
    sscanf(line, "%s %lf %lf %lf", str, c->atoms[i].coor,
                                        c->atoms[i].coor + 1,
                                        c->atoms[i].coor + 2);
    for (int j = 0; j < 3; ++j)
    {
      // TODO two alternatives, one wrapped, one unwrapped coordinates
      c->atoms[i].coor[j] -= floor(c->atoms[i].coor[j]/c->dm[j]) * c->dm[j];
    }
    c->atoms[i].Z = PT_AtomicNumber(sym);
  }

  free(line);
  return 0;
}

int XYZWriteFrame(FILE *f, Crystal *c)
{
  fprintf(f, "%d\ncelldm %10.6f %10.6f %10.6f 90.0 90.0 90.0\n", c->nat, c->dm[0], c->dm[1], c->dm[2]);
  for (int i = 0; i < c->nat; ++i)
  {
    fprintf(f, "%-2s %12.9f %12.9f %12.9f\n", PT_Symbol(c->atoms[i].Z),
                                              c->atoms[i].coor[0],
                                              c->atoms[i].coor[1],
                                              c->atoms[i].coor[2]);
  }
  return 0;
}
