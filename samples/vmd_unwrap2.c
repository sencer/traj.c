// vi: fdm=syntax
#include "xyz.h"
#include "fragments.h"

int vmd_print_out(int *movement, FragmentsInfo *frg, Crystal *c)
{
  FILE *fp;
  int fd;
  char outxyz[40];
  sprintf(outxyz, "/tmp/%s/unw.XXXXXX", getenv("USER"));
  fd = mkstemp(outxyz);
  fp = fdopen(fd, "w");
  int in_flag, *head;
  for (int ax = 0; ax < 3; ++ax)
  {
    for (int dir = 0; dir < 2; ++dir) // iterating twice isn't necessary, two temp files can be used
    {
      in_flag = 0;
      for (int ifrg = 0; ifrg < frg->nfrags; ++ifrg)
      {
        if (movement[ifrg*3+ax]*pow(-1, dir) > 0)
        {
          if(!in_flag)
          {
            in_flag = 1;
            fprintf(fp, "set auto_unwrap_atoms [atomselect top \"index ");
          }
          head = FragmentHead(frg, ifrg);
          for (int iatm = 0; iatm < frg->lfrags[ifrg]; ++iatm)
          {
            fprintf(fp, "%d ", head[iatm]);
          }
        }
      }
      if(in_flag)
      {
        fprintf(fp, "\"]\n");
        fprintf(fp, "$auto_unwrap_atoms moveby \"");
        for (int i = 0; i < ax; ++i)
        {
          fprintf(fp, "0.0 ");
        }
        fprintf(fp, "%14.10f ", pow(-1, dir)*c->dm[ax]);
        for (int i = ax+1; i < 3; ++i)
        {
          fprintf(fp, "0.0 ");
        }
        fprintf(fp, "\"\n");
        fprintf(fp, "$auto_unwrap_atoms delete\n");
      }
    }
  }
  fclose(fp);
  printf("%s", outxyz);
  return 0;
}

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

Crystal *read_file(char *fname)
{
  FILE *f = fopen(fname, "r");
  int nat;
  double dm[3];
  XYZReadHeader(f, &nat, dm);
  Crystal *c =   CrystalInit(nat, dm);
  XYZReadFrame(f, c);
  fclose(f);
  return c;
}

int merge(Crystal *c, BondingInfo *bnd, FragmentsInfo *frg, int x, int *frag_list, int *parent, int *movement)
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

      // calculate the dimensions to move
      for (int k = 0; k < 3; ++k)
      {
        dist = 2 * (c->atoms[iatm1].coor[k] - c->atoms[iatm2].coor[k]);
        if      (dist >  c->dm[k]) { move[k] = -1; }
        else if (dist < -c->dm[k]) { move[k] =  1; }
        movement[ifrg1*3+k] += move[k];
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
          for (int m = 0; m < 3; ++m)
          {
            movement[k*3+m] += move[k];
          }
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
  Crystal *c = read_file(argv[1]);
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
  int parent[frg->nfrags],
  *movement = calloc(3*frg->nfrags, sizeof(int));
  memset(parent, -1, frg->nfrags*sizeof(int));

  for (int i = 0; i < frg->nfrags; ++i)
  {
    merge(c, bnd_p, frg, order[i][0], frag_list, parent, movement);
  }
  free(frag_list);
  BondingDelete(bnd_p);

  vmd_print_out(movement, frg, c);
  FragmentsDelete(frg);
  CrystalDelete(c);
  free(movement);

  return 0;
}
