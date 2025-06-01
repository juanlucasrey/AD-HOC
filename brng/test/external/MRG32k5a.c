/*
   32-bits Random number generator U(0,1): MRG32k5a
   Author: Pierre L'Ecuyer,
   Source: Good Parameter Sets for Combined Multiple Recursive Random
           Number Generators,
           Shorter version in Operations Research,
           47, 1 (1999), 159--164.
   ---------------------------------------------------------
*/
#include "MRG32k5a.h"

#define norm 2.3283163396834613e-10
#define m1   4294949027.0
#define m2   4294934327.0
#define a12     1154721.0
#define a14     1739991.0
#define a15n    1108499.0
#define a21     1776413.0
#define a23      865203.0
#define a25n    1641052.0


/***
The seeds for s10, s11, s12, s13, s14 must be integers in [0, m1 - 1]
   and not all 0. 
The seeds for s20, s21, s22, s23, s24 must be integers in [0, m2 - 1]
   and not all 0. 
***/

#define SEED 12345

static double s10 = SEED, s11 = SEED, s12 = SEED, s13 = SEED, s14 = SEED,
              s20 = SEED, s21 = SEED, s22 = SEED, s23 = SEED, s24 = SEED;


double MRG32k5a (void)
{
   long k;
   double p1, p2;

   /* Component 1 */
   p1 = a12 * s13 - a15n * s10;
   if (p1 > 0.0)
      p1 -= a14 * m1;
   p1 += a14 * s11;
   k = p1 / m1;
   p1 -= k * m1;
   if (p1 < 0.0)
      p1 += m1;
   s10 = s11;
   s11 = s12;
   s12 = s13;
   s13 = s14;
   s14 = p1;

   /* Component 2 */
   p2 = a21 * s24 - a25n * s20;
   if (p2 > 0.0)
      p2 -= a23 * m2;
   p2 += a23 * s22;
   k = p2 / m2;
   p2 -= k * m2;
   if (p2 < 0.0)
      p2 += m2;
   s20 = s21;
   s21 = s22;
   s22 = s23;
   s23 = s24;
   s24 = p2;

   /* Combination */
   if (p1 <= p2)
      return ((p1 - p2 + m1) * norm);
   else
      return ((p1 - p2) * norm);
}