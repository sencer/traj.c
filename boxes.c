#include "boxes.h"

CoarseBox* BoxInit(Crystal *c)
{
  CoarseBox *box = malloc(sizeof(CoarseBox));

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
    box->ngrid[i] = round(c->dm[i]/WIDTH);
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
  int id, comp[3];
  double w[3];

  for (int i = 0; i < 3; ++i)
  {
    w[i] = c->dm[i] / box->ngrid[i];
  }
  for (int i = 0; i < c->nat; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      comp[j] = floor(c->atoms[i].coor[j] / w[j]);
    }
    id = BoxGetIndice(box, comp);
    box->bins[id][box->binsn[id]++]  = i;
  }
  return 0;
}

void GetNeighboringBoxes(CoarseBox *box, int id, int ids[26])
{
    int components1[3],
        components2[3],
        i=0;
    BoxGetComponents(box, id, components1);
    for (int j = 0; j < 27; ++j)
    {
      if(j==13) continue;
      components2[0] = components1[0] + j / 9 - 1;
      components2[1] = components1[1] + (j % 9) / 3 - 1;
      components2[2] = components1[2] + j % 3 - 1;
      ids[i++] = BoxGetIndice(box, components2);
    }
}

int BoxClear(CoarseBox *box)
{
  memset(box->binsn, 0, box->ntot * sizeof(int));
  return 0;
}
