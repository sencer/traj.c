#include "boxes.h"

CoarseBox* BoxInit(Crystal *c, double width)
{
  CoarseBox *box = malloc(sizeof(CoarseBox));
  box->width = width;

  box->bins  = malloc(sizeof(int[MAPB]));
  box->binsn = malloc(sizeof(int));
  box->ntot  = 1;

  BoxUpdate(c, box);

  return box;
}

int BoxUpdate(Crystal *c, CoarseBox *box)
{
  int ntot = 1;

  for (int i = 0; i < 3; ++i)
  {
    box->ngrid[i] = round(c->dm[i]/box->width);
    box->w[i] = c->dm[i] / box->ngrid[i];
    ntot *= box->ngrid[i];
  }
  box->n2d = box->ngrid[1] * box->ngrid[2];

  if (ntot != box->ntot)
  {
    box->ntot = ntot;
    free(box->bins);
    free(box->binsn);
    box->bins  = malloc(ntot * sizeof(int[MAPB]));
    box->binsn = malloc(ntot * sizeof(int));
  }
  BoxClear(box);
  return 0;
}

void BoxDelete(CoarseBox *box)
{
  free(box->bins);
  free(box->binsn);
  free(box);
}

int BoxGetIndice(CoarseBox *box, int comp[3])
{
  return mod(comp[0], box->ngrid[0]) * box->n2d +
         mod(comp[1], box->ngrid[1]) * box->ngrid[2] +
         mod(comp[2], box->ngrid[2]);
}

void BoxGetComponents(CoarseBox *box, int ind, int comp[3])
{
  // TODO Check for negative ind etc.
  comp[0] = ind / box->n2d;
  comp[1] = (ind % box->n2d) / box->ngrid[2];
  comp[2] = ind % box->ngrid[2];
}

int BoxFill(Crystal *c, CoarseBox *box)
{
  int id;

  for (int i = 0; i < c->nat; ++i)
  {
    id = BoxFromCoor(c->atoms[i].coor, box);
    box->bins[id][box->binsn[id]++]  = i;
  }
  return 0;
}

int BoxFromCoor(double coor[3], CoarseBox *box)
{
  int comp[3] = {
    floor(coor[0]/box->w[0]),
    floor(coor[1]/box->w[1]),
    floor(coor[2]/box->w[2])
  };
  return BoxGetIndice(box, comp);
}

void GetNeighboringBoxes(CoarseBox *box, int id, int ids[27])
{
    int components1[3],
        components2[3];
    BoxGetComponents(box, id, components1);
    for (int j = 0; j < 27; ++j)
    {
      components2[0] = components1[0] + j / 9 - 1;
      components2[1] = components1[1] + (j % 9) / 3 - 1;
      components2[2] = components1[2] + j % 3 - 1;
      ids[j] = BoxGetIndice(box, components2);
    }
}

void BoxesOfAtoms(Crystal *c, CoarseBox *box, int *box_of_atom)
{
  // iterate all boxes
  for (int i = 0; i < box->ntot; ++i)
  {
    // and all atoms in the boxes
    for (int j = 0; j < box->binsn[i]; ++j)
    {
      // if the atom is in the target molecule
      box_of_atom[box->bins[i][j]] = i;
    }
  }
}

int BoxClear(CoarseBox *box)
{
  memset(box->binsn, 0, box->ntot * sizeof(int));
  return 0;
}
