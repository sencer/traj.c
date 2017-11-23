#include "xyz.h"
#include "bonding.h"
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

int CheckBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.2||(typ>2 && dist<1.4)||(typ>10 && dist<1.84))?1:0;
}

int main(int argc, char *argv[])
{
  FILE *f = fopen("go.base.xyz", "r");
  int nat, tot_nat=0, oh, epoxy, oh_added = 0, epoxy_added=0, pos_gr=0, pos_fn, oth_atom,
      skip=0, flag=0, nsheets=atoi(argv[3]), *liste;

  double dm[3], tot, dist, zshift=4.8,
         co_ratio = atof(argv[1]),
         epoxy_oh_ratio = atof(argv[2]);

  XYZReadHeader(f, &nat, dm);
  liste = malloc(sizeof(int)*nat);
  for (int i = 0; i < nat; ++i)
  {
    liste[i] = i;
  }
  // Calculate the number of OH and Epoxy groups
  tot = 1.0*nat/co_ratio; // total number of functional groups
  oh = (int) rint(tot/(1+epoxy_oh_ratio));
  epoxy = (int) tot-oh;
  tot_nat = nat+2*oh+epoxy;
  printf("Trying to add %d epoxides and %d OHs to have a C/O ratio %-5.2f\n", epoxy, oh, co_ratio);
  printf("You should have %d atoms if all atoms can be added.\n", nsheets*tot_nat);

  Crystal *c = CrystalInit(tot_nat, dm);
  c->nat = nat;
  XYZReadFrame(f, c);
  fclose(f);
  f = fopen("new_go.xyz", "w");

  // Shuffle the atoms array times.
  srand(getpid()*time(NULL));
  int r, i;

  CoarseBox *box = BoxInit(c, 2.0);
  BoxFill(c, box);

  BondingInfo *bnd = BondingInit(c);
  BondingPopulate(c, box, bnd, CheckBonding, 1);
  BoxDelete(box);

  c->dm[2] = nsheets*zshift;
  c->nat = nsheets*tot_nat;
  XYZWriteHeader(f, c);
  c->nat = nat;
  tot_nat =0;
  for (int sheet = 0; sheet < nsheets; ++sheet)
  {
    for (int j = 0; j < 10*nat; ++j)
    {
      i = j%nat;
      skip = liste[i];
      r = rand()%nat;
      liste[i] = liste[r];
      liste[r] = skip;
    }
    if (argc == 5)
    {
      flag = 1;
    }
    epoxy_added = 0;
    oh_added = 0;
    skip = 0;
    pos_gr = 0;

    while (epoxy_added+oh_added < epoxy+oh && pos_gr<nat+2)
    {
      pos_gr = epoxy_added + oh_added + skip;    // id of current atom on graphene
      pos_fn = nat + epoxy_added + 2 * oh_added; // id of the memory position

      if (flag == 1 && pos_gr == nat)
      {
        printf("Walked over all C atoms, but couldn't add\n");
        printf("all functional groups required. Walking again.\n");
        skip = (epoxy-epoxy_added); // using as a temprorary space!
        oh += skip;
        epoxy = epoxy_added;
        realloc(c->atoms, (2*oh+epoxy+nat)*sizeof(Atom));
        skip = -epoxy_added-oh_added;
        flag++;
        continue;
      }

      // if the current atom is already epoxidized, skip
      if (c->atoms[liste[pos_gr]].id>0)
      {
        skip++;
        continue;
      }

      // try adding epoxy
      r = 0;
      while(r<3&&c->atoms[bnd->bonds[liste[pos_gr]][r]].id>0)
      {
        r++;
      }

      oth_atom = bnd->bonds[liste[pos_gr]][r];
      if (r<3 && epoxy_added<epoxy)
      {
        c->atoms[pos_fn].id = 1;
        c->atoms[pos_fn].Z  = 8;

        dist = c->atoms[oth_atom].coor[0]-c->atoms[liste[pos_gr]].coor[0];
        if(dist>3) { dist -= c->dm[0]; }
        else if (dist<-3) { dist += c->dm[0]; }
        c->atoms[pos_fn].coor[0] = c->atoms[liste[pos_gr]].coor[0] + dist*0.5;

        dist = c->atoms[oth_atom].coor[1]-c->atoms[liste[pos_gr]].coor[1];
        if(dist>3) { dist -= c->dm[1]; }
        else if (dist<-3) { dist += c->dm[1]; }
        c->atoms[pos_fn].coor[1] = c->atoms[liste[pos_gr]].coor[1] + dist*0.5;

        dist = c->atoms[oth_atom].coor[2]-c->atoms[liste[pos_gr]].coor[2];
        if(dist>3) { dist -= c->dm[2]; }
        else if (dist<-3) { dist += c->dm[2]; }
        c->atoms[pos_fn].coor[2] = c->atoms[liste[pos_gr]].coor[2] + dist*0.5;

        r = rand()%2;
        c->atoms[pos_fn].coor[2] += (r?1:-1);

        c->atoms[liste[pos_gr]].id = 1;
        c->atoms[oth_atom].id = 1;

        epoxy_added++;
      }
      else if(oh_added<oh)
      {
        c->atoms[liste[pos_gr]].id = 2;
        r = rand()%2;
        memcpy(c->atoms+pos_fn, c->atoms+liste[pos_gr], sizeof(Atom));
        c->atoms[pos_fn].Z = 8;
        c->atoms[pos_fn].coor[2] += (r?1.4:-1.4);
        memcpy(c->atoms+pos_fn+1, c->atoms+pos_fn, sizeof(Atom));
        c->atoms[1+pos_fn].Z = 1;
        c->atoms[1+pos_fn].coor[2] += (r?0.2:-0.2);
        c->atoms[1+pos_fn].coor[r] += (rand()%2?0.8:-0.8);
        c->atoms[1+pos_fn].coor[(r+1)%2] += (rand()%2?0.5:-0.5);
        oh_added++;
      }
      else
      {
        skip++;
      }
    }
    printf("Added %d epoxides and %d OHs, C/O ratio is %-5.2f\n", epoxy_added, oh_added, nat*1.0/(epoxy_added+oh_added));

    c->nat = nat+epoxy_added+2*oh_added;
    tot_nat += c->nat;
    XYZWriteCoors(f, c);
    c->nat = nat;
    for (int i = 0; i < nat; ++i)
    {
      c->atoms[i].coor[2] += zshift;
      c->atoms[i].id = 0;
    }
  }

  rewind(f);
  fprintf(f, "%d\n", tot_nat);

  free(liste);
  fclose(f);
  BondingDelete(bnd);
  CrystalDelete(c);

  return 0;
}
