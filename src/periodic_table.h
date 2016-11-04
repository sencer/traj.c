#ifndef PERIODIC_TABLE_C_FDEVMGPR
#define PERIODIC_TABLE_C_FDEVMGPR

#include <string.h>
#define PT_HASH_SIZE 483

unsigned char *PT_Symbol(int Z);
unsigned int   PT_AtomicNumber(unsigned char *chr);

#endif /* end of include guard: PERIODIC_TABLE_C_FDEVMGPR */
