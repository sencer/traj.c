#include "xyz.h"
#include "boxes.h"
#include "math.h"
#include "../vendor/hungarian/hungarian.h"

#define MAX_DIST 750 // in picometers
// TODO Needs a larger MAPB in the boxes.h. Use stg like 20.

// for this I don't need to think about pbc. otherwise use CrystDist
int dist(double c1[3], double c2[3])
{
  return 100*sqrt(
    pow((c2[0] - c1[0]), 2) +
    pow((c2[1] - c1[1]), 2) +
    pow((c2[2] - c1[2]), 2)
  );
}

void GetBoxesForAtoms(Crystal *c, CoarseBox *box, int target_nat, int *box_of_atom)
{
  // iterate all boxes
  for (int i = 0; i < box->ntot; ++i)
  {
    // and all atoms in the boxes
    for (int j = 0; j < box->binsn[i]; ++j)
    {
      // if the atom is in the target molecule
      if (box->bins[i][j] < target_nat)
      {
        // note which box it is in.
        box_of_atom[box->bins[i][j]] = i;
      }
    }
  }
}

// TODO: Use this refactoring in Bonding.
void GetCosts(Crystal *c, CoarseBox *box, int atom, int boxid, int target_nat, int costs[target_nat][c->nat-target_nat])
{
  int atom2;

  for (int i = 0; i < box->binsn[boxid]; ++i)
  {
    atom2 = box->bins[boxid][i];
    if ( atom2 >= target_nat)
    {
      if (c->atoms[atom].Z == c->atoms[atom2].Z)
      {
        costs[atom][atom2-target_nat] = dist(c->atoms[atom].coor, c->atoms[atom2].coor);
      }
      else
      {
        costs[atom][atom2-target_nat] = 1000000;
      }
    }
    // For this specific purpose, we are interested into O-O and Ti-Ti
    // distances but not in O-Ti distances; and only when first and second
    // atoms are of "different" molecules'.
    // Normally one should just use atom != atom2
  }
}

int BuildCostMatrix(Crystal *c, CoarseBox *box, int target_nat, int costs[target_nat][c->nat-target_nat])
{
  int *box_of_atom = malloc(sizeof(int) * target_nat), box1, neigh[27];

  GetBoxesForAtoms(c, box, target_nat, box_of_atom);

  for (int atm1 = 0; atm1 < target_nat; ++atm1)
  {
    box1 = box_of_atom[atm1];
    GetNeighboringBoxes(box, box1, neigh);
    for (int j = 0; j < 27; ++j)
    {
      GetCosts(c, box, atm1, neigh[j], target_nat, costs);
    }
  }
  free(box_of_atom);
}

int MapAtoms(Crystal *c, CoarseBox *box, int target_nat, int *atomMap)
{
  int ref_nat = c->nat - target_nat, counter=0;
  // a matrix to hold distances: target_nat rows, reference_nat cols
  int (*costs)[ref_nat] = malloc(sizeof(*costs) * target_nat);

  // fill the matrix with long distances for now.
  for (int i = 0; i < target_nat; ++i)
  {
    for (int j = 0; j < ref_nat; ++j)
    {
      costs[i][j] = 100000;
    }
  }

  BuildCostMatrix(c, box, target_nat, costs);

  hungarian_problem_t p;
  int matrix_size = hungarian_init(&p, target_nat, ref_nat, costs, HUNGARIAN_MODE_MINIMIZE_COST) ;
  /* hungarian_print_costmatrix(&p); */
  hungarian_solve(&p);
  /* hungarian_print_assignment(&p); */

  for (int i = 0; i < target_nat; ++i)
  {
    for (int j = 0; j < ref_nat; ++j)
    {
      if (p.assignment[i][j]==1 && costs[i][j]<MAX_DIST)
      {
        atomMap[i] = target_nat + j;
        atomMap[target_nat + j] = i;
        counter++;
      }
    }
  }

  hungarian_free(&p);
  free(costs);

  return counter;
}

int main(int argc, char *argv[])
{
  if (argc<3) {
    return -1;
  }

  FILE *target = fopen(argv[1], "r");
  FILE *ref = fopen(argv[2], "r");

  int target_nat, ref_nat, nmapped, frame=0;
  double dm[3];


  XYZReadHeader(ref, &ref_nat, dm);
  XYZReadHeader(target, &target_nat, dm);

  // Allocate enough memory to read both structures into one Crystal
  Crystal *c     = CrystalInit(target_nat + ref_nat, dm);
  CoarseBox *box = BoxInit(c, 3.0);

  // this will hold a mapping between the atoms of target and ref molecules
  int *atomMap = malloc(c->nat * sizeof(int));

  // set number of atoms and the c->atoms pointer properly to read ref file
  // the atoms will be read to the end of c->atoms for historical reasons
  c->nat = ref_nat;
  c->atoms = c->atoms + target_nat;
  XYZReadFrame(ref, c);
  fclose(ref);


  // reset pointer position
  c->atoms = (c->atoms - target_nat);

  FILE *out1 = fopen("tar_missing.xyz", "w");
  FILE *out2 = fopen("ref_missing.xyz", "w");

  while(!feof(target))
  {
    printf("\rProcessing frame %d", frame++);
    memset(atomMap, -1, (ref_nat+target_nat) * sizeof(int));
    // read a frame
    c->nat = target_nat;
    XYZReadFrame(target, c);
    c->nat  = target_nat + ref_nat;
    // assign the atoms to boxes
    BoxFill(c, box);
    // Map the atoms in target and ref
    // reset number of atoms and pointer position
    nmapped = MapAtoms(c, box, target_nat, atomMap);
    printf(" %d atoms are mapped.", nmapped);
    fflush(stdout);

    // Print varxyz file.
    fprintf(out1, "%d\ncelldm 100 100 100 90 90 90\n", target_nat - nmapped);
    for (int i = 0; i < target_nat; ++i)
    {
      if(atomMap[i] == -1)
        fprintf(out1, "%2s %10.6f %10.6f %10.6f\n", PT_Symbol(c->atoms[i].Z),
                                                    c->atoms[i].coor[0],
                                                    c->atoms[i].coor[1],
                                                    c->atoms[i].coor[2]);
    }
    fflush(out1);

    // Print varxyz file.
    fprintf(out2, "%d\ncelldm 100 100 100 90 90 90\n", ref_nat - nmapped);
    for (int i = target_nat; i < c->nat; ++i)
    {
      if(atomMap[i] == -1)
        fprintf(out2, "%2s %10.6f %10.6f %10.6f\n", PT_Symbol(c->atoms[i].Z),
                                                    c->atoms[i].coor[0],
                                                    c->atoms[i].coor[1],
                                                    c->atoms[i].coor[2]);
    }
    fflush(out2);

    // read the header information for the next frame
    XYZReadHeader(target, &target_nat, dm);

    // remove the box assignments.
    BoxClear(box);
  }

  fclose(target);
  fclose(out1);
  fclose(out2);

  CrystalDelete(c);
  BoxDelete(box);
  free(atomMap);

  return 0;
}

// vi: fdm=syntax
