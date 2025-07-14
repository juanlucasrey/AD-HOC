/*  Written in 2021 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

#include <stdint.h>
// #include <strings.h>

/* This is a Marsaglia multiply-with-carry generator with period
   approximately 2^255. It is faster than a scrambled linear
   generator, as its only 128-bit operations are a multiplication and sum;
   it is an excellent generator based on congruential arithmetic.

   As all MWC generators, it simulates a multiplicative LCG with prime
   modulus m = 0xfff62cf2ccc0cdaeffffffffffffffffffffffffffffffffffffffffffffffff
   and multiplier given by the inverse of 2^64 modulo m. The modulus has a
   particular form, which creates some theoretical issues, but at this
   size a generator of this kind passes all known statistical tests. For a
   generator of the same type with stronger theoretical guarantees
   consider a Goresky-Klapper generalized multiply-with-carry generator.

   Note that a previous version had a different MWC_A3. Moreover, now
   the result is computed using the current state.
*/

#define MWC_A3 0xfff62cf2ccc0cdaf

/* The state must be initialized so that 0 < c < MWC_A3 - 1.
   For simplicity, we suggest to set c = 1 and x, y, z to a 192-bit seed. */
uint64_t x, y, z, c;

uint64_t inline next() {
	const uint64_t result = z;
	const __uint128_t t = MWC_A3 * (__uint128_t)x + c;
	x = y;
	y = z;
	z = t;
	c = t >> 64;
	return result;
}


/* The following jump functions use a minimal multiprecision library. */

#define MP_SIZE 5
#include "mp.c"

static uint64_t mod[MP_SIZE] = { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, MWC_A3 - 1 };


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int r = cpp_int("0x4b89aa2cd51c37b9f6f8c3fd02ec98fbfe88c291203b225428c3ff11313847eb");
   cpp_int s = ((x + y * b + z * b * b + c * b3) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   z = uint64_t(s >> 128);
   c = uint64_t(s >> 192);
*/

void jump(void) {
	static uint64_t jump[MP_SIZE] = { 0x28c3ff11313847eb, 0xfe88c291203b2254, 0xf6f8c3fd02ec98fb, 0x4b89aa2cd51c37b9 };
	uint64_t state[MP_SIZE] = { x, y, z, c };
	mul(state, jump, mod);
	x = state[0];
	y = state[1];
	z = state[2];
	c = state[3];
}


/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations.
   
   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int r = cpp_int("0xaf5ca22408cdc8306c40ce860e0d702f95382f758ac987764c6e39cf92f77a4");
   cpp_int s = ((x + y * b + z * b * b + c * b3) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   z = uint64_t(s >> 128);
   c = uint64_t(s >> 192);
*/

void long_jump(void) {
	static uint64_t long_jump[MP_SIZE] = { 0x64c6e39cf92f77a4, 0xf95382f758ac9877, 0x6c40ce860e0d702, 0xaf5ca22408cdc83 };
	uint64_t state[MP_SIZE] = { x, y, z, c };
	mul(state, long_jump, mod);
	x = state[0];
	y = state[1];
	z = state[2];
	c = state[3];
}