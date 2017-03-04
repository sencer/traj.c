#include "topo.h"

int *Angles(Crystal *c, BondingInfo *bnd, int *n)
{
  *n = 0;
  for (int i = 0; i < c->nat; ++i)
  {
    for (int j = 0; j < bnd->nbonds[i]-1; ++j)
    {
      for (int k = j+1; k < bnd->nbonds[i]; ++k)
      {
        (*n)++;
      }
    }
  }

  int *angles = (int *)malloc(3 * (*n) * sizeof(int));
  *n = 0;
  for (int i = 0; i < c->nat; ++i)
  {
    for (int j = 0; j < bnd->nbonds[i]-1; ++j)
    {
      for (int k = j+1; k < bnd->nbonds[i]; ++k)
      {
        angles[(*n)*3] = bnd->bonds[i][j];
        angles[(*n)*3+1] = i;
        angles[(*n)*3+2] = bnd->bonds[i][k];
        (*n)++;
      }
    }
  }
  return angles;
}

int *Dihedrals(Crystal *c, BondingInfo *bnd, int *n)
{
  *n = 0;
  int atm1, atm2, atm3, atm4;
  // count the dihedrals
  // for each atom
  for (atm1 = 0; atm1 < c->nat; ++atm1)
  {
    // if the atom has more than one bonds
    if (bnd->nbonds[atm1] > 1)
    {
      // look at each of the bonds (atm1,atm2) is a bond
      for (int j = 0; j < bnd->nbonds[atm1]; ++j)
      {
        atm2 = bnd->bonds[atm1][j];
        // and look at them only once (don't repeat for atm2,atm1)
        // if atm2 also has more than one bonds, then there will be a 
        // dihedral
        if ( atm1 < atm2 && bnd->nbonds[atm2]>1)
        {
          // now for each bond (atm2, atm3)
          for (int k = 0; k < bnd->nbonds[atm2]; ++k)
          {
            atm3 = bnd->bonds[atm2][k];
            // except for the bond (atm2, atm1)
            if (atm3 == atm1) continue;
            // each bond (atm1,atm4) will form a dihedral
            for (int l = 0; l < bnd->nbonds[atm1]; ++l)
            {
              atm4 = bnd->bonds[atm1][l];
              // except for the (atm1,atm2)
              if (atm4 == atm2 || atm4 == atm3) continue;
              (*n)++;
            }
          }
        }
      }
    }
  }

  int *dihedrals = (int *)malloc(4 * (*n) * sizeof(int));
  *n = 0;
  for (atm1 = 0; atm1 < c->nat; ++atm1)
  {
    // if the atom has more than one bonds
    if (bnd->nbonds[atm1] > 1)
    {
      // look at each of the bonds (atm1,atm2) is a bond
      for (int j = 0; j < bnd->nbonds[atm1]; ++j)
      {
        atm2 = bnd->bonds[atm1][j];
        // and look at them only once (don't repeat for atm2,atm1)
        // if atm2 also has more than one bonds, then there will be a 
        // dihedral
        if ( atm1 < atm2 && bnd->nbonds[atm2]>1)
        {
          // now for each bond (atm2, atm3)
          for (int k = 0; k < bnd->nbonds[atm2]; ++k)
          {
            atm3 = bnd->bonds[atm2][k];
            // except for the bond (atm2, atm1)
            if (atm3 == atm1) continue;
            // each bond (atm1,atm4) will form a dihedral
            for (int l = 0; l < bnd->nbonds[atm1]; ++l)
            {
              atm4 = bnd->bonds[atm1][l];
              // except for the (atm1,atm2)
              if (atm4 == atm2 || atm4 == atm3) continue;
              dihedrals[4*(*n)] = atm4;
              dihedrals[4*(*n)+1] = atm1;
              dihedrals[4*(*n)+2] = atm2;
              dihedrals[4*(*n)+3] = atm3;
              (*n)++;
            }
          }
        }
      }
    }
  }
  return dihedrals;
}

int *Impropers(Crystal *c, BondingInfo *bnd,  int *n)
{
  *n = 0;

  for (int i = 0; i < c->nat; ++i)
  {
    if (bnd->nbonds[i] > 2)
    {
      for (int j = 0; j < bnd->nbonds[i]-2; ++j)
      {
        for (int k = j+1; k < bnd->nbonds[i]-1; ++k)
        {
          for (int l = k+1; l < bnd->nbonds[i]; ++l)
          {
            (*n)++;
          }
        }
      }
    }
  }

  int *impropers = (int *)malloc(4 * (*n) * sizeof(int));
  *n = 0;
  for (int i = 0; i < c->nat; ++i)
  {
    if (bnd->nbonds[i] > 2)
    {
      for (int j = 0; j < bnd->nbonds[i]-2; ++j)
      {
        for (int k = j+1; k < bnd->nbonds[i]-1; ++k)
        {
          for (int l = k+1; l < bnd->nbonds[i]; ++l)
          {
            impropers[4*(*n)]   = i;
            impropers[4*(*n)+1] = bnd->bonds[i][j];
            impropers[4*(*n)+2] = bnd->bonds[i][k];
            impropers[4*(*n)+3] = bnd->bonds[i][l];
            (*n)++;
          }
        }
      }
    }
  }
  return impropers;
}

// vi: fdm=syntax
