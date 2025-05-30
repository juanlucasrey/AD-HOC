/*
   64-bits Random number generator U[0,1): lfsr258
   Author: Pierre L'Ecuyer,
   Source: http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
   ---------------------------------------------------------
*/

#include "lfsr258.h"
#include <stdint.h>


/**** VERY IMPORTANT **** :
  The initial seeds y1, y2, y3, y4, y5  MUST be larger than
  1, 511, 4095, 131071 and 8388607 respectively.
****/

#define GERME 123456789123456789

// changed y1 to Y1 because it clashes with some macos definition!!
static uint64_t Y1 = GERME,
    y2 = GERME, y3 = GERME, y4 = GERME, y5 = GERME;


double lfsr258 (void)
{
   uint64_t b;

   b = ((Y1 << 1) ^ Y1) >> 53;
   Y1 = ((Y1 & 18446744073709551614UL) << 10) ^ b;
   b = ((y2 << 24) ^ y2) >> 50;
   y2 = ((y2 & 18446744073709551104UL) << 5) ^ b;
   b = ((y3 << 3) ^ y3) >> 23;
   y3 = ((y3 & 18446744073709547520UL) << 29) ^ b;
   b = ((y4 << 5) ^ y4) >> 24;
   y4 = ((y4 & 18446744073709420544UL) << 23) ^ b;
   b = ((y5 << 3) ^ y5) >> 33;
   y5 = ((y5 & 18446744073701163008UL) << 8) ^ b;
   return (Y1 ^ y2 ^ y3 ^ y4 ^ y5) * 5.421010862427522170037264e-20;
}