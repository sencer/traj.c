#define _GNU_SOURCE
#include "pdb.h"

void trim(char *s)
{
    char *p = s;
    int l = strlen(p);

    while(isspace(p[l - 1])) p[--l] = 0;
    while(* p && isspace(* p)) ++p, --l;

    memmove(s, p, l + 1);
}

// this is not a PDB reader per se. just for personal usage.
int PDBReadHeader(FILE *f, int *nat, double *dm, int first)
{
  char *line = NULL;
  size_t len = 0;
  int read, offset=-1, cur_pos = 0;

  *nat = 0;

  while ( (read = getline(&line, &len, f)) != -1 )
  {
    cur_pos += read;
    if (strncmp("CRYST", line, 5)==0)
    {
      if (offset != -1 ) break;
      sscanf(line, "%*s %lf %lf %lf %*f %*f %*f %*s", dm, dm+1, dm+2);
      if (!first) break;
      offset = cur_pos;
    }
    else if (strncmp("ATOM", line, 4)==0)
    {
      (*nat)++;
    }
  }
  fseek(f, offset, SEEK_SET);
  free(line);

  return 0;
}

int PDBReadFrame(FILE *f, Crystal *c)
{
  int read, counter = 0;
  char *line = NULL, tmp[8];
  size_t len = 0;
  double x, y, z;

  while (counter<c->nat)
  {
    getline(&line, &len, f);
    if (strncmp("ATOM", line, 4)==0)
    {

      memcpy(tmp, &line[30], 8);
      c->atoms[counter].coor[0] = atof(tmp);

      memcpy(tmp, &line[38], 8);
      c->atoms[counter].coor[1] = atof(tmp);

      memcpy(tmp, &line[46], 8);
      c->atoms[counter].coor[2] = atof(tmp);

      memcpy(tmp, &line[76], 2);
      tmp[2] = '\0';
      trim(tmp);
      if (tmp[1]!='\0'&&tmp[1]<91) {
        tmp[1] += 32;
      }
      c->atoms[counter].Z = PT_AtomicNumber(tmp);

      counter++;
    }
  }
  free(line);
  return 0;
}

int PDBWriteFrame(FILE *f, Crystal *c)
{

  fprintf(f, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f-%11s%4d\n", c->dm[0],
      c->dm[1], c->dm[2], 90.0, 90.0, 90.0, "", 0);

  for (int i = 0; i < c->nat; ++i)
  {
    fprintf(f, "ATOM  %5d  %2s      X   1    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
        i+1, PT_Symbol(c->atoms[i].Z), c->atoms[i].coor[0],
        c->atoms[i].coor[1], c->atoms[i].coor[2], 0.0, c->atoms[i].id/100.0,
        PT_Symbol(c->atoms[i].Z));
  }
  fprintf(f, "END\n");
  return 0;
}

