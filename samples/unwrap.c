#include "fragments.h"
#include "lammpstrj.h"
#include "xyz.h"

int checkBonding(double dist, int t1, int t2)
{
  int typ = t1 + t2;
  return (dist<1.0||(typ>2 && dist<1.3)||(typ>10 && dist<2.0))?1:0;
}


int main(int argc, char *argv[])
{
  /* if (argc < 2) { fprintf(stderr, "Reading from the stdin\n"); } */
  int xyz = strcmp(strrchr(argv[1], '.'), ".xyz") != -1;
  FILE *f = (argc>1)?fopen(argv[1], "r"):stdin,
       *f2;
  if (xyz)
    f2 = fopen("unwrap.xyz", "w");
  else
    f2 = fopen("unwrap.lammpstrj", "w");
  int t, nat;
  double dm[3];

  // read the number of atoms and cell dimensions from lammpstrj
  if (xyz)
    XYZReadHeader(f, &nat, dm);
  else
    LMPReadHeader(f, &t, &nat, dm);

  // read the number of atoms and cell dimensions from lammpstrj
  /* LMPReadHeader(f, &t, &nat, dm); */
  // initialize a crystal structure with these dimensions
  Crystal *c = CrystalInit(nat, dm);
  // ...and a coarse graining box
  CoarseBox *box = BoxInit(c, 2.0);
  // ...and a bonding information container
  BondingInfo *bnd = BondingInit(c);
  // ...and a fragmentation information container
  Fragments *frg = FragmentsInit(bnd);

  // now, while we didn't hit the end of file
  while(!feof(f))
  {
    // read a frame
    if (xyz)
      XYZReadFrame(f, c);
    else
      LMPReadFrame(f, c);


    // assign the atoms to boxes
    BoxFill(c, box);
    // populate the bonding list
    BondingPopulate(c, box, bnd, checkBonding);
    // populate the fragments list
    FragmentsPopulate(bnd, frg);
    // Do the work here!
    FragmentsMerge(c, bnd, frg);
    // and write the output

    if (xyz)
      XYZWriteFrame(f2, c);
    else
      LMPWriteFrame(f2, c, t);

    // read the header information for the next frame
    if (xyz)
      XYZReadHeader(f, &nat, dm);
    else
      LMPReadHeader(f, &t, &nat, dm);
    // update cell size
    CrystalSetCell(c, dm);
    // Before moving to the next frame clear the bonding information
    BondingClear(bnd);
    // and the fragmentation information
    FragmentsClear(frg);
    // update the coarse graining box, in case cell dimensions changed
    BoxUpdate(c, box);
  }

  // free the memory
  fclose(f);
  fclose(f2);
  FragmentsClear(frg);
  BondingDelete(bnd);
  BoxDelete(box);
  CrystalDelete(c);

  return 0;
}
