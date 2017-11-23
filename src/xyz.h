#ifndef XYZ_H_CDLMGTJZ
#define XYZ_H_CDLMGTJZ
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "crystal.h"
#include "periodic_table.h"
int XYZReadHeader(FILE *f, int *nat, double *dm);
int XYZReadFrame(FILE *f, Crystal *c);
int XYZWriteFrame(FILE *f, Crystal *c);
int XYZWriteCoors(FILE *f, Crystal *c);
int XYZWriteHeader(FILE *f, Crystal *c);
#endif /* end of include guard: XYZ_H_CDLMGTJZ */
