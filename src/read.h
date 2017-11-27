#ifndef READ_H_XLQV8ZQR
#define READ_H_XLQV8ZQR

#include <stdio.h>
#include "crystal.h"
#include "lammpstrj.h"
#include "xyz.h"

FILE *open_file(char *fname, int *xyz);
Crystal *read_file(FILE *f, int xyz, int *t);
Crystal *read_single_frame(char *fname, int *xyz, int *t);
void *read_next(FILE *f, Crystal *c, int xyz, int *t);

#endif /* end of include guard: READ_H_XLQV8ZQR */
