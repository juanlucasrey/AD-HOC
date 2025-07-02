/*************************************************************************
Copyright:  Francois Panneton and Pierre L'Ecuyer, UniversitÃ© de MontrÃ©al
            Makoto Matsumoto, Hiroshima University

Notice:     This code can be used freely for personal, academic,
            or non-commercial purposes. For commercial purposes,
            please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca

   This code can also be used under the terms of the GNU General Public
   License as published by the Free Software Foundation, either version 3
   of the License, or any later version. See the GPL licence at URL
   http://www.gnu.org/licenses

**************************************************************************/

/* seed is an array of 16 unsigned int's used as seed of the generator */
void InitWELLRNG512a (unsigned int *seed);

double WELLRNG512a (void);