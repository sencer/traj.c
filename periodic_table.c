#include "periodic_table.h"

// This is to implement a simple hash table for { element_sym: atomic_number }
// returns a unique integer for each element symbol
static unsigned int hash(unsigned char *chr)
{
  if (strlen((char*)chr) == 1)
  {
    return chr[0] + 353;
  }
  else
  {
    return (chr[0]-66) + 20 * ( chr[1]-97 );
  }
}

unsigned char *PT_Symbol(int Z)
{
 return ELEMENTS[Z-1];
}

unsigned int  PT_AtomicNumber(unsigned char *chr)
{
  const int hash_table[] = {
    56, 20, 0, 0, 0, 31, 0, 0, 0, 0, 57, 0, 11, 0, 91, 0, 88, 0, 73, 0, 0, 0,
    105, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41, 0, 82, 0, 37, 51, 65, 89, 0, 0, 0, 70,
    0, 0, 0, 0, 0, 0, 0, 115, 0, 0, 0, 0, 0, 21, 43, 0, 0, 48, 0, 0, 0, 64, 0,
    0, 0, 0, 0, 101, 60, 0, 46, 0, 0, 0, 0, 0, 4, 58, 0, 0, 26, 32, 2, 0, 0,
    0, 0, 0, 10, 0, 0, 0, 75, 34, 52, 0, 0, 98, 54, 0, 0, 0, 72, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 104, 0, 0, 47, 0, 0, 0, 0, 0, 0, 80, 0, 0, 0, 0, 12, 0,
    118, 0, 0, 111, 106, 0, 0, 107, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 113, 0,
    0, 0, 45, 0, 90, 0, 83, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 28, 0, 0, 0, 0,
    14, 22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 97,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13, 0, 17, 0, 0,
    114, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 81, 95, 0, 96, 0, 0, 100, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 61, 0, 0, 62, 69, 0, 0, 112, 0, 0, 0, 0, 0, 49, 0,
    0, 0, 25, 0, 0, 0, 0, 86, 50, 0, 0, 0, 27, 0, 0, 30, 0, 67, 0, 0, 0, 0,
    42, 102, 0, 84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 93, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    18, 35, 24, 0, 68, 87, 0, 0, 77, 0, 36, 103, 0, 0, 0, 59, 0, 0, 38, 0, 33,
    0, 55, 110, 99, 40, 0, 108, 0, 0, 0, 0, 0, 0, 76, 0, 0, 0, 0, 117, 85, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 109, 0, 0, 78, 0, 0, 0, 0, 79, 0, 29, 0, 63,
    0, 0, 0, 0, 0, 0, 71, 0, 0, 0, 94, 0, 44, 0, 0, 5, 6, 0, 0, 9, 0, 1, 53,
    0, 19, 0, 116, 7, 8, 15, 0, 0, 16, 0, 92, 23, 74, 0, 39, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 65
  };

  return hash_table[hash(chr)];
}

// vi: fdm=syntax
