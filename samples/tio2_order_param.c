// vi: fdm=syntax
#include "boxes.h"
#include "read.h"

#define PI 3.14159265
#define RAD2DEG 180/PI

#define UNKNOWN     0
#define ANATASE     1
#define RUTILE      2
#define TWIN    3
#define TETRAHEDRAL 4

#define THR 10

#define TI_1ST_NEIGH 10.5625 //11.2225
#define TI_2ND_NEIGH 15.21
#define O_NEIGH      5.29

int Common(int *xlist, int xlen, int *ylist, int ylen)
{
  // receive to lists ( & their lengths), and return 1 if they have
  // a common element
  for (int i = 0; i < xlen; ++i)
  {
    for (int j = 0; j < ylen; ++j)
    {
      if (xlist[i] == ylist[j])
      {
        return 1;
      }
    }
  }
  return 0;
}

int TriangleType(int x, int y, int n1[][5], int l1[], int n2[][10], int l2[])
{
  // Search for triangles as described in doi:10.1021/jp301228x
  // for anatase, rutile and anatase {112} twins
  //
  // Will receive the first and second nearest neighbors of two Ti atoms
  // (x and y) which themselves are first degree nearest neighbors

  if (Common(n1[x], l1[x], n2[y], l2[y]) || Common(n2[x], l2[x], n1[y], l1[y]))
  {
    // if the first nearest neighbors of one Ti atom and the second nearest
    // neighbor of the other is the same atom, triangle type 1 exists that is a
    // triangle with 2 edges = nearest neighbor distance and the third edge is
    // second nearest neighbor distance
    return 1;
  }
  else if (Common(n2[x], l2[x], n2[y], l2[y]))
  {
    // or else if the two Ti atoms has a common second-nearest-neighbor
    // triangle is type 2; that is it has 2 edges with 2nd nearest neigbor
    // distance, and one with nearest neighbor distance
    return 2;
  }
  // no such triangle exists
  return 0;
}

int double_cmp(const void *a, const void *b)
{
  // compare two doubles
  const double *fa = (const double *)a;
  const double *fb = (const double *)b;
  return (int)(1000*(*fa-*fb));
}

void VecDiff(double u[3], double v[3], double diff[3])
{
  // difference of two 3D vectors
  for (int i = 0; i < 3; ++i)
  {
    diff[i] = u[i] - v[i];
  }
}

double VecLen2(double u[3])
{
  // length square of a 3D vector
  return pow(u[0], 2) + pow(u[1], 2) + pow(u[2], 2);
}

double Dist2(double u[3], double v[3])
{
  // squared euclidian distance between two 3D vectors
  double diff[3];
  VecDiff(u, v, diff);
  return VecLen2(diff);
}

double VecDot(double u[3], double v[3])
{
  // dot product of two 3D vectors
  double ret = 0;

  for (int i = 0; i < 3; ++i)
  {
    ret += u[i] * v[i];
  }

  return ret;
}

double VecAngle(double u[3], double v[3])
{
  // angle between two 3D vectors, in degrees
  double lu = sqrt(VecLen2(u)),
         lv = sqrt(VecLen2(v));
  return RAD2DEG * acos(VecDot(u, v) / (lu * lv));
}

int main(int argc, char *argv[])
{
  int t, xyz, nat;
  FILE *f = open_file(argv[1], &xyz);
  Crystal *c = read_file(f, xyz, &t);
  CoarseBox *box = Box(c, 4.00); // enough to ensure second-nearest neigbors not missed

  // construct a map and inverse map of Ti atoms, so that
  // c->atom[ map[i] ] will give the i-th Ti atom and
  // invmap[map[i]] will give i.
  // this should save some memory

  int nti = 0, // number of Ti atoms
      map[c->nat/2], // TiO2 has 1/3 Ti, but leave room for some reduction
      invmap[c->nat];

  // count Ti atoms, and build the maps
  for (int i = 0; i < c->nat; ++i)
  {
    if (c->atoms[i].Z == 22)
    {
      map[nti] = i;
      invmap[i] = nti;
      nti++;
    }
  }

  // construct tables of nearest / second nearest Ti atoms, and nearest O atoms
  // for each Ti

  int ti_nearest[nti][5], num_ti_nearest[nti],
  ti_2nearest[nti][10], num_ti_2nearest[nti],
  o_nearest[nti][8], num_o_nearest[nti],
  ii, // ii = map[i]
  jj; // jj is the id of the other atom

  // will iterate through boxes, some aux data
  int current_box, other_box, neighbor_box[27], num_neighbor_box,
      box_list[c->nat];

  // output will written to
  FILE *out = fopen("annotated.xyz", "w");

  // do it for all the frames in the input
  int counter = 1, eof = 0;
  do
  {
    // populate the box_list, so box_list[i] will give the id of the box 
    // of the i-th atom
    BoxesOfAtoms(c, box, box_list);
    // reset nearest neighbor data
    memset(num_ti_nearest,  0, nti*sizeof(int));
    memset(num_ti_2nearest, 0, nti*sizeof(int));
    memset(num_o_nearest,   0, nti*sizeof(int));

    // loop through all the Ti atoms in the frame
    for (int i = 0; i < nti; ++i)
    {
      ii = map[i];
      c->atoms[ii].id = 0;
      current_box = box_list[ii];

      num_neighbor_box = GetNeighboringBoxes(box, current_box, neighbor_box, 0);

      // iterate through the neighboring boxes (which include the box itself)
      for (int j = 0; j < num_neighbor_box; ++j)
      {
        other_box = neighbor_box[j];
        // and all the atoms in the neighboring box
        for (int k = 0; k < box->binsn[other_box]; ++k)
        {
          jj = box->bins[other_box][k];
          if (c->atoms[jj].Z == 22)
          {
            if ( jj <= ii ) continue; // don't count Ti - Ti twice

            if (Dist2(c->atoms[ii].coor, c->atoms[jj].coor) < TI_1ST_NEIGH)
            {
              ti_nearest[i][num_ti_nearest[i]++] = jj;
              ti_nearest[invmap[jj]][num_ti_nearest[invmap[jj]]++] = ii;
            }
            else if (Dist2(c->atoms[ii].coor, c->atoms[jj].coor) < TI_2ND_NEIGH)
            {
              ti_2nearest[i][num_ti_2nearest[i]++] = jj;
              ti_2nearest[invmap[jj]][num_ti_2nearest[invmap[jj]]++] = ii;
            }
          }
          else if ( c->atoms[jj].Z == 8 &&
              Dist2(c->atoms[ii].coor, c->atoms[jj].coor) < O_NEIGH)
          {
            o_nearest[i][num_o_nearest[i]++] = jj;
          }
        }
      }
    }

    // decide on the type of Ti atoms
    int triangle, all_zero, all_one, all_two;
    double angle, u[3], v[3], w[3], x[3], angles[6];
    for (int i = 0; i < nti; ++i)
    {
      ii = map[i];
      if (num_ti_nearest[i] > 4)
      {
        c->atoms[ii].id = UNKNOWN;
        continue;
      }
      all_zero = 1;
      all_one = 1;
      all_two = 1;
      // i is one corner of the triangle
      for (int j = 0; j < num_ti_nearest[i]; ++j)
      {
        jj = ti_nearest[i][j]; // jj is the second corner
        triangle = TriangleType(i, invmap[jj],
            ti_nearest, num_ti_nearest, ti_2nearest, num_ti_2nearest
            );
        if ( triangle == 1)
        {
          all_zero = 0;
          all_two = 0;
        }
        else if (triangle == 2)
        {
          all_one = 0;
          all_zero = 0;
        }
      }

      /* now decide on Ti type according to 10.1021/jp301228x
       * 0 -> unknown
       * 1 -> anatase / brookite
       * 2 -> rutile
       * 3 -> anatase twin
       * 4 -> tetrahedral
       */
      if (all_zero)
      {
        c->atoms[ii].id = UNKNOWN;
      }
      if (all_one)
      {
        c->atoms[ii].id = ANATASE;
      }
      else if (num_ti_nearest[i] == 3 && !all_two && !all_zero)
      {
        c->atoms[ii].id = TWIN;
      }
      else if (num_ti_nearest[i] == 2 && all_two)
      {
        VecDiff(c->atoms[ti_nearest[i][0]].coor, c->atoms[ii].coor, u);
        VecDiff(c->atoms[ti_nearest[i][1]].coor, c->atoms[ii].coor, v);
        angle = VecAngle(u, v);
        if ( angle > 155 )
        {
          c->atoms[ii].id = RUTILE;
        }
        else
        {
          c->atoms[ii].id = TWIN;
        }
      }
      else
      {
        c->atoms[ii].id = UNKNOWN;
      }
      if (num_o_nearest[i] == 4)
      {
        VecDiff(c->atoms[ii].coor, c->atoms[o_nearest[i][0]].coor, u);
        VecDiff(c->atoms[ii].coor, c->atoms[o_nearest[i][1]].coor, v);
        VecDiff(c->atoms[ii].coor, c->atoms[o_nearest[i][2]].coor, w);
        VecDiff(c->atoms[ii].coor, c->atoms[o_nearest[i][3]].coor, x);
        angles[0] = VecAngle(u, v);
        angles[1] = VecAngle(u, w);
        angles[2] = VecAngle(u, x);
        angles[3] = VecAngle(v, w);
        angles[4] = VecAngle(v, x);
        angles[5] = VecAngle(w, x);
        qsort(angles, 6, sizeof(double), double_cmp);

        if (angles[0] + angles[1] + angles[2] + angles[3] > 395)
        {
          c->atoms[ii].id = TETRAHEDRAL;
        }
      }
    }
    // do a second pass over UNKNOWN Ti with one
    // nearest neighbor and print
    fprintf(out, "%d\ncelldm 100 100 100\n", nti);
    for (int i = 0; i < nti; ++i)
    {
      ii = map[i];
      if ( c->atoms[ii].id == UNKNOWN && num_ti_nearest[i] == 1)
      {
        c->atoms[ii].id = c->atoms[ti_nearest[i][0]].id;
      }
      fprintf(out, "Ti %14.10f %14.10f %14.10f %d\n", c->atoms[ii].coor[0],
          c->atoms[ii].coor[1], c->atoms[ii].coor[2], c->atoms[ii].id); }

    printf("%d\n", counter);
    for (int i = 0; i < 50 && eof != EOF; ++i)
    {
      counter++;
      read_next(f, c, xyz, &t);
      eof = fgetc(f);
      fseek(f, -1, SEEK_CUR);
    }
    BoxUpdate(c, box);
    BoxFill(c, box);
  }
  while(eof != EOF);

  fclose(out);
  fclose(f);

  BoxDelete(box);
  CrystalDelete(c);
  return 0;
}
