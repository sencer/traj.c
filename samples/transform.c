// vi: fdm=syntax
#include <unistd.h>
#include "read.h"
#include "hungarian.h"

#define CENTER 75.0
#define MAX_BETA 48.0
#define DO_TRANSFORM 1

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


int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    printf("np_rmsd.c traj.xyz\n");
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

  Crystal *ref = read_file(ftraj, xyz, &t);
  Transform(ref);

  do
  {
    for (int i = 0; i < 1; ++i)
    {
      eof = fgetc(ftraj);
      if (eof == EOF) break;
      fseek(ftraj, -1, SEEK_CUR);
      read_next(ftraj, ref, xyz, &t);
      Transform(ref);
      counter++;
      printf("\rProcessing frame %d", counter);
      fflush(stdout);
    }

    LMPWriteFrame(fout, ref, t);

  }
  while (eof!=EOF);

  printf("\n");

  CrystalDelete(ref);
  fclose(ftraj);
  fclose(fout);

  return 0;
}
