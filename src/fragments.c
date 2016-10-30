#include "fragments.h"

Fragments *FragmentsInit(BondingInfo *bnd)
{
  Fragments *frg = malloc(sizeof(Fragments));

  frg->frags  = malloc(bnd->nat * sizeof(int));
  frg->lfrags = malloc(bnd->nat * sizeof(int));
  frg->nfrags = 0;

  BondingClear(bnd);

  return frg;
}

void FragmentsDelete(Fragments *frg)
{
  free(frg->frags);
  free(frg->lfrags);
  free(frg);
}

int FragmentsClear(Fragments *frg)
{
  frg->nfrags = 0;
  return 0;
}

int FragmentsPopulate(BondingInfo *bnd, Fragments *frg)
{
  int visited[bnd->nat], // a list to mark the visited atoms
      lvisit[bnd->nat],  // a list of not-yet-visited atoms in the fragment
      pvisit = 0,        // current position in the lvisit array
      pfrag = 0,         // current position in the frg->frags array
      cur,               // "current atom"
      ldiff[MBPA],       // a list store difference of two bonding lists
      lend;              // a place in memory to store length of ldiff

  memset(visited, 0, bnd->nat * sizeof(int));
  // we will iterate through each atom, to get to which atoms they are bound
  for (int i = 0; i < bnd->nat; ++i)
  {
    if (!visited[i]) // skip if we already visited this atom
    {
      frg->lfrags[frg->nfrags] = 1 + bnd->bondsn[i];
      frg->frags[pfrag++] = i; // add the atom i to current fragment, and move
      visited[i] = 1;

      // now find each atom that are bound to i and add them to
      for (int j = 0; j < bnd->bondsn[i]; ++j)
      {
        lvisit[pvisit++] = bnd->bonds[i][j];    // to-visit list
        frg->frags[pfrag++] = bnd->bonds[i][j]; // add to current fragment
        visited[bnd->bonds[i][j]] = 1;          // and mark visited.
      }

      // now, while I have more atoms to visit in my lvisit list
      while(pvisit)
      {
        // get the last atom from the lvisit list, and move head backwards
        cur = lvisit[--pvisit];
        // then get the list of atoms the "cur" atom is bonded to but
        // not already included in the fragments list to ldiff
        lend = 0;
        for (int j = 0; j < bnd->bondsn[cur]; ++j)
        {
          if(!visited[bnd->bonds[cur][j]])
          {
            ldiff[lend++] = bnd->bonds[cur][j];
          }
        }

        for (int j = 0; j < lend; ++j)
        {
          lvisit[pvisit++] = ldiff[j];    // put those atoms to-visit list
          frg->frags[pfrag++] = ldiff[j]; // and to current fragment
          visited[ldiff[j]] = 1;          // and mark as visited
        }
        frg->lfrags[frg->nfrags] += lend;
      }
      // save the length of fragment
      frg->nfrags++;
    }
  }
  return 0;
}

int FragmentsMerge(Crystal *c, BondingInfo *bnd, Fragments *frg)
{
  double center[3], *c1, *c2, dist;
  int visited[c->nat], cur = 0, oth = 0, counter = 0;
  memset(visited, -1, c->nat*sizeof(int));

  // Iterate through each fragment
  for (int i = 0; i < frg->nfrags; ++i)
  {
    // set fragment center to 0
    memset(center, 0, 3*sizeof(double));

    // iterate through each atom in the fragment
    for (int j = 0; j < frg->lfrags[i]; ++j)
    {
      // frags hold the atoms list in order of fragments they belong
      cur = frg->frags[counter]; // the atom we are looking at

      visited[cur] = 0; // it was -1
      c1 = c->atoms[cur].coor;
      // add the coordinates to the center vector.
      for(int k = 0; k < 3; ++k)
      {
        center[k] += c1[k];
      }
      // iterate through all atoms cur is bonded to 
      for (int k = 0; k < bnd->bondsn[cur]; ++k)
      {
        oth = bnd->bonds[cur][k]; // oth is the atom cur is bonded to
        if (visited[oth] == -1)   // skip if we already did this atom.
        {
          c2 = c->atoms[oth].coor;
          for (int m = 0; m < 3; ++m)
          {
            dist = 2 * (c2[m] - c1[m]); // if distance is more than half of cell
            if      (dist >  c->dm[m]) { c2[m] -= c->dm[m]; }
            else if (dist < -c->dm[m]) { c2[m] += c->dm[m]; }
          }
        }
      }
      counter++;
    }
    for (int j = 0; j < 3; ++j)
    {
      center[j] = floorf(center[j]/frg->lfrags[i]/c->dm[j]);
      if (center[j] == 0) { continue; }
      center[j] *= c->dm[j];
      for (int k = counter-1; k > counter-1-frg->lfrags[i]; --k)
      {
        c->atoms[frg->frags[k]].coor[j] -= center[j];
      }
    }
  }
  return 0;
}
