// vi: fdm=syntax
#include <unistd.h>
#include "read.h"
#include "bonding.h"

#define CENTER 50.0

double beta(double coor[3])
{
  double center[3] = {CENTER, CENTER, CENTER};

  return sqrt(pow(coor[0]-center[0], 2) +
              pow(coor[1]-center[1], 2) +
              pow(coor[2]-center[2], 2));
}

int check_bonding(double len, int t1, int t2)
{
  int t = t1 + t2;
  return ((len < 1.7 && t > 10) ||
          (len < 2.3 && t > 22));
  /* return ( len < 1.2 || */
  /*         (len < 1.4 && t > 2) || */
  /*         (len < 1.7 && t > 10) || */
  /*         (len < 2.3 && t > 22)); */
}


int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    printf("coordination traj.{xyz,lammpstrj}\n");
    return 1;
  }

  if (access(argv[1], F_OK) == -1)
  {
    printf("%s does not exist.\n", argv[1]);
    return 1;
  }

  // start reading trajectory
  int t, xyz, nat, eof, counter = 1;

  FILE *ftraj = open_file(argv[1], &xyz),
       *fout  = fopen("out.xyz", "w");

  Crystal *c = read_file(ftraj, xyz, &t);
  CoarseBox *box = Box(c, 2.4);
  BondingInfo *bnd = Bonds(c, box, check_bonding, 0);
  Atom at;

  nat = 0;
  for (int i = 0; i < c->nat; ++i)
  {
    at = c->atoms[i];
    if (at.Z > 1)
    {
      nat++;
    }
  }

  do
  {
    eof = fgetc(ftraj);
    fprintf(fout, "%d\ncelldm %f %f %f\n", nat, 2*CENTER, 2*CENTER, 2*CENTER);
    for (int i = 0; i < c->nat; ++i)
    {
      at = c->atoms[i];
      if (at.Z > 1)
      {
        fprintf(fout, "%-2s %8.4f %8.4f %8.4f %6.2f %d\n", PT_Symbol(at.Z),
            at.coor[0], at.coor[1], at.coor[2], beta(at.coor), bnd->nbonds[i]);
      }
    }
    read_next(ftraj, c, xyz, &t);
    BoxClear(box);
    BoxFill(c, box);
    BondingClear(bnd);
    BondingPopulate(c, box, bnd, check_bonding, 0);
  }
  while (eof!=EOF);

  printf("\n");

  CrystalDelete(c);
  BoxDelete(box);
  fclose(ftraj);
  fclose(fout);

  return 0;
}
