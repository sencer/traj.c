#include "read.h"

FILE *open_file(char *fname, int *xyz)
{
  FILE *f;

  if ((int)fname[0] == '-')
  {
    fprintf(stderr, "Reading from the stdin\n");
    f = stdin;
    *xyz = !strncmp(fname+1, "xyz", 3);
  }
  else
  {
    f = fopen(fname, "r");
    *xyz = strcmp(strrchr(fname, '.'), ".xyz") != -1;
  }
  return f;
}

Crystal *read_file(FILE *f, int xyz, int *t)
{
  Crystal *c;
  int nat;
  double dm[3];
  if (xyz)
  {
    XYZReadHeader(f, &nat, dm);
    c = CrystalInit(nat, dm);
    XYZReadFrame(f, c);
  }
  else
  {
    LMPReadHeader(f, t, &nat, dm);
    c = CrystalInit(nat, dm);
    LMPReadFrame(f, c);
  }
  return c;
}

Crystal *read_single_frame(char *fname, int *xyz, int *t)
{
  FILE *f = open_file(fname, xyz);
  Crystal *c = read_file(f, *xyz, t);
  fclose(f);
  return c;
}
