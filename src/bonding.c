#include "bonding.h"

BondingInfo *Bonds(Crystal *c, CoarseBox *box, int (*cb)(double, int, int), int periodic)
{
  BondingInfo *bnd = BondingInit(c);
  BondingPopulate(c, box, bnd, cb, periodic);
  return bnd;
}

BondingInfo *BondingInit(Crystal *c)
{
  BondingInfo *bnd = malloc(sizeof(BondingInfo));

  bnd->nat = c->nat;
  bnd->bonds  = malloc(bnd->nat * sizeof(int[MBPA]));
  bnd->nbonds = calloc(bnd->nat, sizeof(int));

  return bnd;
}

void BondingDelete(BondingInfo *bnd)
{
  free(bnd->bonds);
  free(bnd->nbonds);
  free(bnd);
}

int BondingDoBonding(Crystal *c, BondingInfo *bnd, int i, int j, int (*cb)(double, int, int))
{
  double dist = CrystalDist(c, c->atoms[i].coor, c->atoms[j].coor);

  if(cb(dist, c->atoms[i].Z, c->atoms[j].Z))
  {
    bnd->bonds[i][bnd->nbonds[i]++] = j;
    bnd->bonds[j][bnd->nbonds[j]++] = i;
  }
  return 0;
}

int BondingPopulate(Crystal *c, CoarseBox *box, BondingInfo *bnd, int (*cb)(double, int, int), int periodic)
{
  int icomp[3], jcomp[3], jbox, atm1, atm2, neigh[27], numneigh;

  // for each box
  for (int ibox = 0; ibox < box->ntot; ++ibox)
  {
    // if there are atoms in the box
    if (box->binsn[ibox])
    {
      // cross-check bonding between the atoms in box ibox
      for (int i = 0; i < box->binsn[ibox]-1; ++i)
      {
        atm1 = box->bins[ibox][i];
        for (int j = i+1; j < box->binsn[ibox]; ++j)
        {
          atm2 = box->bins[ibox][j];
          BondingDoBonding(c, bnd, atm1, atm2, cb);
        }
      }
      // now check the neighboring boxes - only the 13 of them in (+) direction
      numneigh = GetNeighboringBoxes(box, ibox, neigh, periodic);
      BoxGetComponents(box, ibox, icomp);
      for (int j = 0; j < numneigh; ++j)
      {
        jbox = neigh[j];
        if (ibox >= jbox) { continue; } // don't do boxes twice
        // for each atom in box jbox check bonding
        for (int k = 0; k < box->binsn[jbox]; ++k)
        {
          atm1 = box->bins[jbox][k];
          for (int i = 0; i < box->binsn[ibox]; ++i)
          {
            atm2 = box->bins[ibox][i];
            BondingDoBonding(c, bnd, atm1, atm2, cb);
          }
        }
      }
    }
  }
  return 0;
}

int BondingClear(BondingInfo *bnd)
{
  memset(bnd->nbonds, 0, bnd->nat * sizeof(int));
  return 0;
}

int BondingPrint(BondingInfo *bnd)
{
  for (int i = 0; i < bnd->nat; ++i)
  {
    printf("%-4d: ", i);
    for (int j = 0; j < bnd->nbonds[i]; ++j)
    {
      printf("%d ", bnd->bonds[i][j]);
    }
    printf("\n");
  }
  return 0;
}

void BondingReadFile(char fname[], BondingInfo *bnd)
{
  FILE *f = fopen(fname, "r");
  char *line = NULL, *field;
  size_t len = 0;
  int read, offset=-1, cur_pos = 0, atom = 0;

  while ( (read = getline(&line, &len, f)) != -1 )
  {
    field = strtok(line, " ");

    while (field != NULL)
    {
      bnd->bonds[atom][bnd->nbonds[atom]] = atoi(field);
      bnd->nbonds[atom]++;
      field = strtok (NULL, " ");
    }
    atom++;
  }
  free(field);
  free(line);
}

int BondingBFSWalk(BondingInfo *bnd, int node, int *frag, int *checked)
{
  int current, other, lfrag = 0, nqueued =0,
      *queue = calloc(bnd->nat, sizeof(int));

  queue[nqueued++]  = node;
  checked[node] = 1;

  while (nqueued>0)
  {
    current = queue[--nqueued];
    frag[lfrag++] = current;

    for (int j = 0; j < bnd->nbonds[current]; ++j)
    {
      other = bnd->bonds[current][j];
      if ( !checked[other] )
      {
        checked[other] = 1;
        queue[nqueued++] = other;
      }
    }
  }
  free(queue);
  return lfrag;
}
