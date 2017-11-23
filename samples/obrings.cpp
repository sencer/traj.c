extern "C" {
  #include <lammpstrj.h>
  #include <fragments.h>
}
#include <openbabel/mol.h>
/* #include <openbabel/obconversion.h> */
#define SHEET 50

using namespace std;

/*
 * Reads a trajectory to analyze the functional groups on FGSs, and the small
 * molecules formed from these functional groups.
 */

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0 ||
         (typ>2 && dist<1.4) ||
         (typ>10 && dist<1.9) ) ? 1:0;
         /* (t1==6 && t2==6 && dist<2.1))?1:0; */
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); }
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin;
  int t, nat, pos, id, nrings[4], size;
  double dm[3];

  LMPReadHeader(f, &t, &nat, dm);
  Crystal *c = CrystalInit(nat, dm);
  CoarseBox *box = BoxInit(c, 2.0);
  BondingInfo *bnd = BondingInit(c);
  FragmentsInfo *frg = FragmentsInit(bnd);

  int mapping[600] , atom_id;

  OpenBabel::OBMol mol;
  OpenBabel::OBAtom atom;
  vector<OpenBabel::OBRing*> rings;
  vector<OpenBabel::OBRing*>::iterator it;
  /* vector<int>::iterator jt; */
  /* vector<OpenBabel::OBRing*> *rings; */

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    // read a frame
    LMPReadFrame(f, c);
    BoxFill(c, box);
    BondingPopulate(c, box, bnd, checkBonding, 1);
    FragmentsPopulate(bnd, frg);

    pos = 0;
    memset(nrings, 0, 4*sizeof(int));
    for (int i = 0; i < frg->nfrags; ++i)
    {
      if(frg->lfrags[i] > SHEET)
      {
        atom_id = 1;
        memset(mapping,-1,600);

        for (int j = 0; j < frg->lfrags[i]; ++j)
        {
          id = frg->frags[pos++];
          mapping[id] = atom_id++;
          /* if (id<471) */
          /* { */
          atom.SetAtomicNum(c->atoms[id].Z);
          atom.SetVector(c->atoms[id].coor[0], c->atoms[id].coor[1],
              c->atoms[id].coor[2]);
          mol.AddAtom(atom);
          /* } */
        }
        pos -= frg->lfrags[i];
        for (int j = 0; j < frg->lfrags[i]; ++j)
        {
          id = frg->frags[pos++];
          for (int k = 0; k < bnd->nbonds[id]; ++k)
          {
            if (id < bnd->bonds[id][k])
            {
              mol.AddBond(mapping[id], mapping[bnd->bonds[id][k]], 1);
            }
          }
        }
        /* printf("%d %d\n", atom_id, mol.NumBonds()); */
        rings = mol.GetSSSR();

        for (it = rings.begin(); it != rings.end(); ++it)
        {
          size = (*it)->Size();
          if ( size > 4 && size < 9)
          {
            nrings[size-5]++;
            /* if (size == 5) */
            /* { */
            /* for ( jt = (*it)->_path.begin(); jt != (*it)->_path.end(); ++jt) */
            /* { */
            /* printf("%d ", *jt); */
            /* } */
            /* printf("\n"); */
            /* } */
          }
        }
        /* free(mapping); */
      }
    }

    printf("%10d %6d %6d %6d %6d\n", t, nrings[0], nrings[1], nrings[2],
        nrings[3]);
    LMPReadHeader(f, &t, &nat, dm);
    if (nat != c->nat)
    {
      CrystalDelete(c);
      BondingDelete(bnd);
      FragmentsDelete(frg);
      c = CrystalInit(nat, dm);
      bnd = BondingInit(c);
      frg = FragmentsInit(bnd);
    }
    else
    {
      CrystalSetCell(c, dm);
      BondingClear(bnd);
      FragmentsClear(frg);
    }
    fflush(stdout);
    mol.Clear();
    BoxUpdate(c, box);
  }

  // free the memory
  fclose(f);
  FragmentsDelete(frg);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  /* std::ofstream ofs("test.cml"); */
  /* OpenBabel::OBConversion ob(NULL, &ofs); */
  /* ob.SetOutFormat("CML"); */
  /* ob.Write(&mol); */
  /* ob.CloseOutFile(); */

  return 0;
}
