#define _GNU_SOURCE
#include "lammpstrj.h"

int XYZReadHeader(FILE *f, int *nat, double *dm)
{
  char *line = NULL;
  size_t len = 0;
  // TODO Check the section headers being skipped now

  getline(&line, &len, f);
  sscanf(line, "%d", nat);
  getline(&line, &len, f);
  sscanf(line, "%*s %lf %lf %lf %*f %*f %*f", dm, dm+1, dm+2);

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
      c->atoms[i].coor[j] -= floor(c->atoms[i].coor[j]/c->dm[j]) * c->dm[j];
    }
    if     (strcmp(str, "H" ) == 0) { c->atoms[i].Z = 1;  }
    else if(strcmp(str, "C" ) == 0) { c->atoms[i].Z = 6;  }
    else if(strcmp(str, "N" ) == 0) { c->atoms[i].Z = 7;  }
    else if(strcmp(str, "O" ) == 0) { c->atoms[i].Z = 8;  }
    else if(strcmp(str, "Ti") == 0) { c->atoms[i].Z = 22; }
  }

  return 0;
}

// TODO two alternatives, one wrapped, one unwrapped coordinates
