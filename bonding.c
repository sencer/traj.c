#include "bonding.h"

BondingInfo *BondingInit(Crystal *c)
{
  BondingInfo *bnd = malloc(sizeof(BondingInfo));

  bnd->nat = c->nat;
  bnd->bonds  = malloc(bnd->nat * sizeof(int[MBPA]));
  bnd->bondsn = malloc(bnd->nat * sizeof(int));
  bnd->frags  = malloc(bnd->nat * sizeof(int));
  bnd->lfrags = malloc(bnd->nat * sizeof(int));
  bnd->nfrags = 0;

  BondingClear(bnd);

  return bnd;
}

void BondingDelete(BondingInfo *bnd)
{
  free(bnd->frags);
  free(bnd->lfrags);
  free(bnd->bonds);
  free(bnd->bondsn);
  free(bnd);
}

int BondingDoBonding(Crystal *c, BondingInfo *bnd, int i, int j, int (*cb)(double, int, int))
{
  double dist = CrystalDist(c, i, j);
  if(cb(dist, c->atoms[i].Z, c->atoms[j].Z))
  {
    bnd->bonds[i][bnd->bondsn[i]++] = j;
    bnd->bonds[j][bnd->bondsn[j]++] = i;
  }
  return 0;
}

int BondingPopulate(Crystal *c, CoarseBox *box, BondingInfo *bnd, int (*cb)(double, int, int))
{
  int icomp[3], jcomp[3], jbox, atm1, atm2;

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
      BoxGetComponents(box, ibox, icomp);
      for (int j = 14; j < 27; ++j)
      {
        jcomp[0] = icomp[0] + j / 9 - 1;
        jcomp[1] = icomp[1] + (j % 9) / 3 - 1;
        jcomp[2] = icomp[2] + j % 3 - 1;
        jbox = BoxGetIndice(box, jcomp);
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
  bnd->nfrags = 0;
  memset(bnd->bondsn, 0, bnd->nat * sizeof(int));
  return 0;
}

int BondingPrint(BondingInfo *bnd)
{
  for (int i = 0; i < bnd->nat; ++i)
  {
    printf("%-4d: ", i);
    for (int j = 0; j < bnd->bondsn[i]; ++j)
    {
      printf("%d ", bnd->bonds[i][j]);
    }
    printf("\n");
  }
  return 0;
}

int BondingFragments(BondingInfo *bnd)
{
  int visited[bnd->nat], // a list to mark the visited atoms
      lvisit[bnd->nat],  // a list of not-yet-visited atoms in the fragment
      pvisit = 0,        // current position in the lvisit array
      pfrag = 0,         // current position in the bnd->frags array
      cur,               // "current atom"
      ldiff[MBPA],       // a list store difference of two bonding lists
      lend;              // a place in memory to store length of ldiff

  memset(visited, 0, bnd->nat * sizeof(int));
  // we will iterate through each atom, to get to which atoms they are bound
  for (int i = 0; i < bnd->nat; ++i)
  {
    // skip if we already visited this atom
    if (!visited[i])
    {
      bnd->lfrags[bnd->nfrags] = 1 + bnd->bondsn[i];
      bnd->frags[pfrag++] = i; // add the atom i to current fragment, and move

      // now find each atom that are bound to i and add them to
      for (int j = 0; j < bnd->bondsn[i]; ++j)
      {
        lvisit[pvisit++] = bnd->bonds[i][j];    // to-visit list
        bnd->frags[pfrag++] = bnd->bonds[i][j]; // add to current fragment
        visited[bnd->bonds[i][j]] = 1;          // and mark visited.
      }

      // now, while I have more atoms to visit in my lvisit list
      while(pvisit)
      {
        // get the last atom from the lvisit list, and move head backwards
        cur = lvisit[--pvisit];
        // then get the list of atoms the "cur" atom is bonded to but
        // not already included in the fragments list to ldiff
        listDiff(bnd->bonds[cur], bnd->bondsn[cur],
                 bnd->frags,      pfrag,
                 ldiff,           &lend);

        for (int j = 0; j < lend; ++j)
        {
          lvisit[pvisit++] = ldiff[j];    // put those atoms to-visit list
          bnd->frags[pfrag++] = ldiff[j]; // and to current fragment
          visited[ldiff[j]] = 1;          // and mark as visited
        }
        bnd->lfrags[bnd->nfrags] += lend;

      }
      // save the length of fragment
      bnd->nfrags++;
    }
  }
  return 0;
}

int BondingMergeFragments(Crystal *c, BondingInfo *bnd)
{
  double center[3], *c1, *c2, dist;
  int visited[c->nat], cur = 0, oth = 0, counter = 0;
  memset(visited, -1, c->nat*sizeof(int));

  for (int i = 0; i < bnd->nfrags; ++i)
  {
    memset(center, 0, 3*sizeof(int));
    for (int j = 0; j < bnd->lfrags[i]; ++j)
    {
      cur = bnd->frags[counter];
      visited[cur] = 0;
      c1 = c->atoms[cur].coor;
      for(int k = 0; k < 3; ++k)
      {
        center[k] += c1[k];
      }
      for (int k = 0; k < bnd->bondsn[cur]; ++k)
      {
        oth = bnd->bonds[cur][k];
        if (visited[oth] == -1)
        {
          c2 = c->atoms[oth].coor;
          for (int m = 0; m < 3; ++m)
          {
            dist = 2 * (c2[m] - c1[m]);
            if      (dist >  c->dm[m]) { c2[m] -= c->dm[m]; }
            else if (dist < -c->dm[m]) { c2[m] += c->dm[m]; }
          }
        }
      }
      counter++;
    }
    for (int j = 0; j < 3; ++j)
    {
      center[j] = ((int)((center[j]/bnd->lfrags[i])/c->dm[j]))*c->dm[j];
      for (int k = 0; k < bnd->lfrags[i]; ++k)
      {
        c->atoms[counter-1-k].coor[j] -= center[j];
      }
    }
  }
  return 0;
}
