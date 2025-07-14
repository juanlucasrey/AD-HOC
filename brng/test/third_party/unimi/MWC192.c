/*  Written in 2023 by Sebastiano Vigna (vigna@acm.org)

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
#include <strings.h>

/* This is a Marsaglia multiply-with-carry generator with period
   approximately 2^191. It is faster than a scrambled linear
   generator, as its only 128-bit operations are a multiplication and sum;
   it is an excellent generator based on congruential arithmetic.

   As all MWC generators, it simulates a multiplicative LCG with prime
   modulus m = 0xffa04e67b3c95d85ffffffffffffffffffffffffffffffff
   and multiplier given by the inverse of 2^64 modulo m. The modulus has a
   particular form, which creates some theoretical issues, but at this
   size a generator of this kind passes all known statistical tests. For a
   generator of the same type with stronger theoretical guarantees
   consider a Goresky-Klapper generalized multiply-with-carry generator.
*/

#define MWC_A2 0xffa04e67b3c95d86

/* The state must be initialized so that 0 < c < MWC_A2 - 1.
   For simplicity, we suggest to set c = 1 and x, y to a 128-bit seed. */
uint64_t x, y, c;

uint64_t inline next() {
	const uint64_t result = y;
	const __uint128_t t = MWC_A2 * (__uint128_t)x + c;
	x = y;
	y = t;
	c = t >> 64;
	return result;
}


/* The following jump functions use a minimal multiprecision library. */

#define MP_SIZE 4
#include "mp.c"

static uint64_t mod[MP_SIZE] = { 0xffffffffffffffff, 0xffffffffffffffff, MWC_A2 - 1 };

/* This is the jump function for the generator. It is equivalent
   to 2^96 calls to next(); it can be used to generate 2^96
   non-overlapping subsequences for parallel computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b2 = b * b;
   cpp_int m = MWC_A2 * b2 - 1;
   cpp_int r = cpp_int("0xdc2be36e4bd21a2afc217e3b9edf985d94fb8d87c7c6437");
   cpp_int s = ((x + y * b + c * b2) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   c = uint64_t(s >> 128);
*/

void jump(void) {
   static uint64_t jump[MP_SIZE] = { 0xd94fb8d87c7c6437, 0xafc217e3b9edf985, 0xdc2be36e4bd21a2 };
   uint64_t state[MP_SIZE] = { x, y, c };
	mul(state, jump, mod);
	x = state[0];
	y = state[1];
	c = state[2];
}

/* This is the long-jump function for the generator. It is equivalent to
   2^144 calls to next(); it can be used to generate 2^96 starting points,
   from each of which jump() will generate 2^96 non-overlapping
   subsequences for parallel distributed computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b2 = b * b;
   cpp_int m = MWC_A2 * b2 - 1;
   cpp_int r = cpp_int("0x3c6528aaead6bbddec956c3909137b2dd0e7cedd16a0758e");
   cpp_int s = ((x + y * b + c * b2) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   c = uint64_t(s >> 128);
*/

void long_jump(void) {
   static uint64_t long_jump[MP_SIZE] = { 0xd0e7cedd16a0758e, 0xec956c3909137b2d, 0x3c6528aaead6bbdd };
   uint64_t state[MP_SIZE] = { x, y, c };
	mul(state, long_jump, mod);
	x = state[0];
	y = state[1];
	c = state[2];
}