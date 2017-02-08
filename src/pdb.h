#ifndef PDB_H_FDNLALSA
#define PDB_H_FDNLALSA
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "crystal.h"
#include "periodic_table.h"
int PDBReadFrame(FILE *f, Crystal *c);
int PDBReadHeader(FILE *f, int *nat, double *dm, int first);
int PDBWriteFrame(FILE *f, Crystal *c);
#endif /* end of include guard: PDB_H_FDNLALSA */
