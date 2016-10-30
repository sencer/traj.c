#ifndef LAMMPSTRJ_H_WG6ZNYS3
#define LAMMPSTRJ_H_WG6ZNYS3
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "crystal.h"
#include "periodic_table.h"
int LMPReadHeader(FILE *f, int *t, int *nat, double *dm);
int LMPReadFrame(FILE *f, Crystal *c);
int LMPWriteFrame(FILE *f, Crystal *c, int t);
#endif /* end of include guard: LAMMPSTRJ_H_WG6ZNYS3 */
