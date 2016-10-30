#include "util.h"

int isInList(int elem, int *list, int llen)
{
  for (int i = 0; i < llen; ++i)
  {
    if (list[i] == elem)
      return 1;
  }
  return 0;
}

int listDiff(int *l1, int len1, int *l2, int len2, int *ld, int *lend)
{
  // difference of l1 with len1 elements, and l2 with len2 elements
  // are written into ld with lend elements
  // enough memory for ld should be allocated in the caller
  *lend = 0;
  for (int i = 0; i < len1; ++i)
  {
    if (!isInList(l1[i], l2, len2)) ld[(*lend)++] = l1[i];
  }
  return 0;
}

int mod(int x, int y)
{
  int m = x % y;
  return m < 0 ? m+y : m;
}
