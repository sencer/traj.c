#include "xyz.h"
#include "bonding.h"

/* Find the rings in a structure, and output the ring sizes */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}

int getIndex(int needle, int *haystack)
{
  int i = 0;
  while(haystack[i] > -1) {
    if (haystack[i] == needle)
    {
      return i;
    }
    i++;
  }
  return -1;
}

void  printstack(int *stack, int len)
{
  fprintf(stderr,"%d: %d", len, stack[0]);
  for (int i = 1; i < len; ++i)
  {
    fprintf(stderr,"-%d", stack[i]);
  }
  fprintf(stderr,"\n");
}

int main(int argc, char *argv[])
{
  FILE *f = fopen(argv[1], "r");
  int nat;
  double dm[3];

  XYZReadHeader(f, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  XYZReadFrame(f, c);
  CoarseBox *box = BoxInit(c, 2.0);
  BoxFill(c, box);
  BondingInfo *bnd = BondingInit(c);
  BondingPopulate(c, box, bnd, checkBonding, 1);

  // do a depth first search
  int visited[c->nat], stack[c->nat], dpos[c->nat], slen=0, cur=-1, index=-1;
  memset(visited, 0, c->nat*sizeof(int));
  memset(stack, -1, c->nat*sizeof(int));
  memset(dpos, 0, c->nat*sizeof(int));
  // visit an atom
  for (int i = 0; i < c->nat; ++i)
  {
    if (visited[i])
    {
      continue;
    }
    visited[i] = 1;
    stack[slen++] = i;
        /* fprintf(stderr, "NEW FOR I\n"); */

    while (slen>0){
      /* fprintf(stderr, "NEW WHILE\n"); */
      for (int j = dpos[slen-1]; j < bnd->nbonds[stack[slen-1]]; ++j)
      {
        /* fprintf(stderr, "NEW FOR J\n"); */
        /* printstack(stack, slen); */
        cur = bnd->bonds[stack[slen-1]][j];
        /* fprintf(stderr, "Looking at bond %d of atom %d; %d.\n", j, stack[slen-1], cur); */
        if (cur == stack[slen-2]) {
          /* fprintf(stderr, "Not using: cur: %d, grandpa: %d\n", cur, stack[slen-2]); */
          continue;
        }
        index = getIndex(cur, stack);
        if ( index > -1 && slen-index-1 > 1)
        {
          printf("Found a ring of size %d ::  %d", slen-index, stack[index]);
          for (int k = 0; k < slen-index-1; ++k)
          {
            printf(" %d", stack[index+k+1]);
          }
          printf("\n");
        }
        else if (visited[cur])
        {
          visited[cur] = 1;
          dpos[slen-1] = j+1;
          j = -1;
          stack[slen++] = cur;
        }
      }
      stack[--slen] = -1;
    }
  }

  // free the memory
  fclose(f);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
