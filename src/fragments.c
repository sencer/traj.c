// vi: fdm=syntax
#include "fragments.h"

FragmentsInfo *Fragments(BondingInfo *bnd)
{
  FragmentsInfo *frg = FragmentsInit(bnd);
  FragmentsPopulate(bnd, frg);
  return frg;
}

FragmentsInfo *FragmentsInit(BondingInfo *bnd)
{
  FragmentsInfo *frg = malloc(sizeof(FragmentsInfo));

  frg->nat = bnd->nat;

  frg->frags  = malloc(bnd->nat * sizeof(int));
  frg->lfrags = malloc(bnd->nat * sizeof(int));
  frg->nfrags = 0;

  return frg;
}

void FragmentsDelete(FragmentsInfo *frg)
{
  free(frg->frags);
  free(frg->lfrags);
  free(frg);
}

int FragmentsClear(FragmentsInfo *frg)
{
  frg->nfrags = 0;
  return 0;
}

int FragmentsPopulate(BondingInfo *bnd, FragmentsInfo *frg)
{
  int *checked = calloc(bnd->nat, sizeof(int)),
      nfrag = 0, lfrag = 0,
      pos = 0;

  for (int i = 0; i < bnd->nat; ++i)
  {
    if(checked[i]) continue;
    lfrag = BondingBFSWalk(bnd, i, frg->frags+pos, checked);
    pos += lfrag;
    frg->lfrags[nfrag] = lfrag;
    nfrag++;
  }
  free(checked);
  frg->nfrags = nfrag;
  return 0;
}

int FragmentsMerge(Crystal *c, BondingInfo *bnd, FragmentsInfo *frg)
{

  double center[3], *c1, *c2, dist; //, cell_center[3];
  int visited[c->nat], cur = 0, oth = 0, counter = 0;
  memset(visited, -1, c->nat*sizeof(int));

  /* for (int i = 0; i < 3; ++i) */
  /* { */
  /* cell_center[i] = c->dm[i] * 0.5; */
  /* } */

  // Iterate through each fragment
  for (int i = 0; i < frg->nfrags; ++i)
  {
    // set fragment center to 0
    /* memset(center, 0, 3*sizeof(double)); */

    // sort atoms with respect to their distance to the cell center
    // first allocate atom_dist struct
    /* atom_dist = malloc(frg->lfrags[i]*sizeof(AtomDist)); */
    /* for (int j = 0; j < frg->lfrags[i]; ++j) */
    /* { */
    /*   // calculate the distance to cell center */
    /*   cur = frg->frags[counter++]; // the atom we are looking at */
    /*   atom_dist[j].dist = CrystalDist(c, c->atoms[cur].coor, cell_center); */
    /*   atom_dist[j].id = cur; */
    /* } */

    /* qsort(atom_dist, frg->lfrags[i], sizeof(AtomDist), cmpfunc); */

    // iterate through each atom in the fragment
    for (int j = 0; j < frg->lfrags[i]; ++j)
    {
      // frags hold the atoms list in id of fragments they belong
      /* cur = atom_dist[j].id; // the atom we are looking at */
      cur = frg->frags[counter++]; // the atom we are looking at

      visited[cur] = 0; // it was -1
      c1 = c->atoms[cur].coor;
      // add the coordinates to the center vector.
      for(int k = 0; k < 3; ++k)
      {
        center[k] += c1[k];
      }
      // iterate through all atoms cur is bonded to 
      for (int k = 0; k < bnd->nbonds[cur]; ++k)
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

    /* free(atom_dist); */
  }
  return 0;
}

int* FragmentHead(FragmentsInfo *frg, int n)
{
  int l = 0;
  for (int i = 0; i < n; ++i)
  {
    l += frg->lfrags[i];
  }
  return frg->frags + l;
}

int* FragmentList(FragmentsInfo *frg)
{
  int *list = malloc(frg->nat * sizeof(int)),
      pos = 0;

  for (int i = 0; i < frg->nfrags; ++i)
  {
    for (int j = 0; j < frg->lfrags[i]; ++j)
    {
      list[frg->frags[pos++]] = i;
    }
  }
  return list;
}
