#include "bonding.h"

BondingInfo *BondingInit(Crystal *c)
{
  BondingInfo *bnd = malloc(sizeof(BondingInfo));

  bnd->nat = c->nat;
  bnd->bonds  = malloc(bnd->nat * sizeof(int[MBPA]));
  bnd->bondsn = malloc(bnd->nat * sizeof(int));
  bnd->frags  = malloc(2*bnd->nat*sizeof(int));
  bnd->nfrags = 0;

  BondingClear(bnd);

  return bnd;
}

void BondingDelete(BondingInfo *bnd)
{
  free(bnd->frags);
  free(bnd->bonds);
  free(bnd->bondsn);
  free(bnd);
}

int BondingPopulate(Crystal *c, CoarseBox *box, BondingInfo *bnd)
{
  int list[13*MAPB],
      icomp[3],
      jcomp[3],
      lpos,
      jbox,
      atm1,
      atm2;

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
          if(checkBonding(c, atm1, atm2))
          {
            bnd->bonds[atm1][bnd->bondsn[atm1]++] = atm2;
            bnd->bonds[atm2][bnd->bondsn[atm2]++] = atm1;
          }
        }
      }
      // now create a list of atoms in the neighboring boxes
      lpos = 0;
      BoxGetComponents(box, ibox, icomp);
      for (int j = 14; j < 27; ++j) // 14-27, just looking forward.
      {
        jcomp[0] = icomp[0] + j / 9 - 1;
        jcomp[1] = icomp[1] + (j % 9) / 3 - 1;
        jcomp[2] = icomp[2] + j % 3 - 1;
        jbox = BoxGetIndice(box, jcomp);
        if(box->binsn[jbox])
        {
          memcpy(list+lpos, box->bins[jbox], box->binsn[jbox]*sizeof(int));
          lpos += box->binsn[jbox];
        }
      }
      // and check bonding between the atoms in ibox and the neighboring boxes
      for (int i = 0; i < box->binsn[ibox]; ++i)
      {
        atm1 = box->bins[ibox][i];
        for (int j = 0; j < lpos; ++j)
        {
          atm2 = list[j];
          if (checkBonding(c, atm1, atm2))
          {
            bnd->bonds[atm1][bnd->bondsn[atm1]++] = atm2;
            bnd->bonds[atm2][bnd->bondsn[atm2]++] = atm1;
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

int checkBonding(Crystal *c, int i, int j)
{
  double dist = CrystalDist(c, i, j);
  int typ = c->atoms[i].Z + c->atoms[j].Z;
  // TODO this line is very system specific
  // Perhaps a function pointer might be used for BondingPopulate?
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<1.8))?1:0;
}

int BondingFragments(BondingInfo *bnd)
{
  int lvisit[bnd->nat], // a list of not-yet-visited atoms in the fragment
      pvisit = 0,     // current position in the lvisit array
      pfrag = 0,      // current position in the bnd->frags array
      cur,            // "current atom"
      ldiff[MBPA],    // a list store difference of two bonding lists
      lend;           // a place in memory to store length of ldiff

  // we will iterate through each atom, to get to which atoms they are bound
  for (int i = 0; i < bnd->nat; ++i)
  {
    // skip if we already visited this atom
    if (!isInList(i, bnd->frags, pfrag))
    {
      bnd->frags[pfrag++] = i; // add the atom i to current fragment, and move

      // append the atoms bound to i to lvisit array
      memcpy(lvisit+pvisit, bnd->bonds[i], bnd->bondsn[i] * sizeof(int));
      pvisit += bnd->bondsn[i]; // and move the head

      // and also to the bnd->frags
      memcpy(bnd->frags+pfrag, bnd->bonds[i], bnd->bondsn[i] * sizeof(int));
      pfrag += bnd->bondsn[i]; // and move the head

      // now, while I have more atoms to visit in my lvisit list
      while(pvisit)
      {
        // get the last atom from the lvisit list, and move head backwards
        cur = lvisit[--pvisit];
        // then get the list of atoms the "cur" atom is bonded to but
        // not already included in the fragments list to ldiff
        listDiff(bnd->bonds[cur], bnd->bondsn[cur],
            bnd->frags,      pfrag-1,
            ldiff,           &lend);

        // now append these atoms in ldiff to lvisit
        memcpy(lvisit+pvisit, ldiff, lend * sizeof(int));
        pvisit += lend; // and move the head

        // and to the current fragment
        memcpy(bnd->frags+pfrag, ldiff, lend * sizeof(int));
        pfrag += lend; // and move the head
      }
      // if I don’t have anymore atoms to visit, the current fragment is over
      // mark it with a -1
      bnd->frags[pfrag++] = -1;
      bnd->nfrags++;
    }
  }
  return 0;
}
