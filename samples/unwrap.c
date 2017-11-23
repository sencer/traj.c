// vi: fdm=syntax
#include "read.h"
#include "fragments.h"

int check_bonding(double len, int t1, int t2)
{
  int t = t1 + t2;
  return (
      len < 1.2 ||
     (len < 1.4 && t > 2) ||
     (len < 1.7 && t > 10) ||
     (len < 2.3 && t > 22)
  );
}

int frag_sort_len(const void *p, const void *q)
{
  int *x = (int *)p;
  int *y = (int *)q;
  return (x[1] > y[1]);
}

int merge(Crystal *c, BondingInfo *bnd, FragmentsInfo *frg, int x, int *frag_list, int *parent)
{
  int *head = FragmentHead(frg, x), iatm1, iatm2, ifrg1, ifrg2, new_par, par2,
      move[3] = {0, 0, 0};
  double dist;

  for (int i = 0; i < frg->lfrags[x]; ++i)
  {
    iatm1 = head[i];
    ifrg1 = frag_list[iatm1];
    for (int j = 0; j < bnd->nbonds[iatm1]; ++j)
    {
      iatm2 = bnd->bonds[iatm1][j];
      ifrg2 = frag_list[iatm2];
      par2  = parent[ifrg2];

      // skip if the two atoms are in the same fragment already, or
      // they have been linked
      if ( ifrg1 == ifrg2 || ifrg1 == par2 ) continue;

      // otherwise merging into ifrg2, or to its parent if it has
      new_par = ((par2 == -1) ? ifrg2 : par2);
      /* printf("Merging %d into %d through %d.\n", ifrg1, new_par, ifrg2); */

      // calculate the dimensions to move
      for (int k = 0; k < 3; ++k)
      {
        dist = 2 * (c->atoms[iatm1].coor[k] - c->atoms[iatm2].coor[k]);
        if      (dist >  c->dm[k]) { move[k] = -1; }
        else if (dist < -c->dm[k]) { move[k] =  1; }
      }

      // update ifrg1
      parent[ifrg1] = new_par;
      for (int l = 0; l < frg->lfrags[x]; ++l)
      {
        for (int k = 0; k < 3; ++k)
        {
          c->atoms[head[l]].coor[k] += move[k]*c->dm[k];
        }
      }

      // update its children
      for (int k = 0; k < frg->nfrags; ++k)
      {
        if (parent[k] == ifrg1)
        {
          parent[k] = new_par;
          head = FragmentHead(frg, k);
          for (int l = 0; l < frg->lfrags[k]; ++l)
          {
            for (int m = 0; m < 3; ++m)
            {
              c->atoms[head[l]].coor[m] += move[m]*c->dm[m];
            }
          }
        }
      }
      return 0;
    }
  }
  return 1;
}

int main(int argc, char *argv[])
{
  int t, xyz;
  Crystal *c = read_single_frame(argv[1], &xyz, &t);

  CoarseBox *box = Box(c, 2.3);
  BondingInfo *bnd_p = Bonds(c, box, check_bonding, 1);
  BondingInfo *bnd_n = Bonds(c, box, check_bonding, 0);
  BoxDelete(box);
  FragmentsInfo *frg = Fragments(bnd_n);
  BondingDelete(bnd_n);

  int *frag_list = FragmentList(frg);

  // sort the fragments by length
  int order[frg->nfrags][2];
  for (int i = 0; i < frg->nfrags; ++i)
  {
    order[i][0] = i;
    order[i][1] = frg->lfrags[i];
  }
  qsort(order, frg->nfrags, sizeof(order[0]), frag_sort_len);

  // now start from the smallest fragment, merge them into a "parent"
  int parent[frg->nfrags];
  memset(parent, -1, frg->nfrags*sizeof(int));

  for (int i = 0; i < frg->nfrags; ++i)
  {
    merge(c, bnd_p, frg, order[i][0], frag_list, parent);
  }
  free(frag_list);
  BondingDelete(bnd_p);
  FragmentsDelete(frg);

  FILE  *fo;
  if(xyz)
  {
    fo = fopen("unwrap.xyz", "w");
    XYZWriteFrame(fo, c);
  }
  else
  {
    fo = fopen("unwrap.lammpstrj", "w");
    LMPWriteFrame(fo, c, t);
  }

  fclose(fo);

  CrystalDelete(c);
  return 0;
}
