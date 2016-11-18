#include "boxes.h"
#include "xyz.h"
#include "../vendor/hungarian/hungarian.h"

void MapAtoms(Crystal *ref, CoarseBox *ref_box, Crystal *tar, CoarseBox *tar_box)
{
  int (*costs)[ref->nat] = malloc( ref->nat*sizeof(*costs) ),
      r_box_list[ref->nat], r_box_id, tar_atom, neigh[27];

  // get an array, that keeps which atom (array index) is in which box
  // (array value): list[atom_index] = box_index_the_atom_is_in
  for (int i = 0; i < ref_box->ntot; ++i)
  {
    for (int j = 0; j < ref_box->binsn[i]; ++j)
    {
      r_box_list[ref_box->bins[i][j]] = i;
    }
  }

  // create an initial cost matrix for the atoms in the reference frame
  // vs current frame. initial costs are high.
  for (int i = 0; i < ref->nat; ++i)
  {
    for (int j = 0; j < ref->nat; ++j)
    {
      costs[i][j] = 100000;
    }
  }

  // for each atom in the reference frame
  for (int ref_atom = 0; ref_atom < ref->nat; ++ref_atom)
  {
    // get which box is the atom in
    r_box_id = r_box_list[ref_atom];
    // get the ids of the boxes that are neighboring to that box
    GetNeighboringBoxes(ref_box, r_box_id, neigh);

    // iterate these boxes (actually the same boxes for the current frame,
    // ---not ref_box, but tar_box--)
    for (int neigh_box = 0; neigh_box < 27; ++neigh_box)
    {
      // to get all the atoms in each of them
      for (int i = 0; i < tar_box->binsn[neigh_box]; ++i)
      {
        tar_atom = tar_box->bins[neigh_box][i];
        if(ref->atoms[ref_atom].Z == tar->atoms[tar_atom].Z)
        {
          // and write the distance between these atoms and the reference atom.
          costs[ref_atom][tar_atom] = (int) ( 100 * sqrt (
            pow( ref->atoms[ref_atom].coor[0] - tar->atoms[tar_atom].coor[0], 2)+
            pow( ref->atoms[ref_atom].coor[1] - tar->atoms[tar_atom].coor[1], 2)+
            pow( ref->atoms[ref_atom].coor[2] - tar->atoms[tar_atom].coor[2], 2)
          ));
        }
        else
        {
          // if the atoms are of different elements, there is no way the two will
          // be a match. so write a very high cost.
          costs[ref_atom][tar_atom] = 1000000;
        }
      }
    }
  }

  // solve the hungarian problem for bipartite matching
  hungarian_problem_t p;
  int matrix_size = hungarian_init(&p, ref->nat, ref->nat, costs, HUNGARIAN_MODE_MINIMIZE_COST) ;
  hungarian_solve(&p);

  for (int i = 0; i < ref->nat; ++i)
  {
    for (int j = 0; j < tar->nat; ++j)
    {
      if ( p.assignment[i][j] == 1)
      {
        // save the mapping between the reference and target atoms to the
        // id place holder of the reference atoms
        ref->atoms[i].id = j;
//      tar->atoms[j].id = i;
      }
    }
  }

  hungarian_free(&p);
  free(costs);
}

int main(int argc, char *argv[])
{

  if (argc < 2)
  {
    printf("Usage: ./sort_frames traj.xyz\n");
    return -1;
  }
  FILE *out = fopen("sorted.xyz", "w");


  FILE *traj = fopen(argv[1], "r");
  int nat;
  double dm[3];


  // read the first frame as reference frame.
  XYZReadHeader(traj, &nat, dm);
  Crystal *ref = CrystalInit(nat, dm);
  XYZReadFrame(traj, ref);
  CoarseBox *ref_box = BoxInit(ref, 3.0);
  BoxFill(ref, ref_box);

  // and write it as is to the output
  XYZWriteFrame(out, ref);

  // get ready for the following frames.
  Crystal *frm = CrystalInit(nat, dm);
  CoarseBox *frm_box = BoxInit(frm, 3.0);
  int tar_atom;

  while(!feof(traj))
  {
    // this is just to skip the header.
    XYZReadHeader(traj, &nat, dm);

    //read the frame and bin the elements
    XYZReadFrame(traj, frm);
    BoxFill(frm, frm_box);

    // solve the problem.
    MapAtoms(ref, ref_box, frm, frm_box);

    // write the sorted frame
    fprintf(out, "%d\ncelldm %10.6f %10.6f %10.6f 90.0 90.0 90.0\n", frm->nat,
        frm->dm[0], frm->dm[1], frm->dm[2]);

    for (int i = 0; i < frm->nat; ++i)
    {
      tar_atom = ref->atoms[i].id; // this is the sorting data.
      fprintf(out, " %-4s%12.9f %12.9f %12.9f\n",
          PT_Symbol(frm->atoms[tar_atom].Z),
          frm->atoms[tar_atom].coor[0],
          frm->atoms[tar_atom].coor[1],
          frm->atoms[tar_atom].coor[2]);
    }

    fgetc(traj);  // to check if we hit EOF
  }

  fclose(traj);
  fclose(out);
  CrystalDelete(ref);
  CrystalDelete(frm);
  BoxDelete(ref_box);
  BoxDelete(frm_box);
  return 0;
}

// vi: fdm=syntax
