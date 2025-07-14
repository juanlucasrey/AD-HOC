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
// #include <strings.h>

/* This is a Marsaglia multiply-with-carry generator with period
   approximately 2^127. It is close in speed to a scrambled linear
   generator, as its only 128-bit operations are a multiplication and sum;
   it is an excellent generator based on congruential arithmetic.

   As all MWC generators, it simulates a multiplicative LCG with
   prime modulus m = 0xffebb71d94fcdaf8ffffffffffffffff and
   multiplier given by the inverse of 2^64 modulo m. The modulus has a
   particular form, which creates some theoretical issues, but at this
   size a generator of this kind passes all known statistical tests,
   except for very large size birthday-spacings tests in three dimensions
   (as predicted from the theory). These failures are unlikely to have
   any impact in practice.

   For a generator of the same type with stronger theoretical guarantees,
   consider a Goresky-Klapper generalized multiply-with-carry generator.

   Alternatively, applying a xorshift to the result (i.e., returning
   x ^ (x << 32) instead of x) will eliminate the birthday-spacings
   failures.

   Note that a previous version had a different MWC_A1. Moreover, now
   the result is computed using the current state.
*/

#define MWC_A1 0xffebb71d94fcdaf9

/* The state must be initialized so that 0 < c < MWC_A1 - 1.
   For simplicity, we suggest to set c = 1 and x to a 64-bit seed. */
uint64_t x, c;

uint64_t inline next() {
	const uint64_t result = x; // Or, result = x ^ (x << 32) (see above)
	const __uint128_t t = MWC_A1 * (__uint128_t)x + c;
	x = t;
	c = t >> 64;
	return result;
}


/* The following jump functions use a minimal multiprecision library. */

#define MP_SIZE 3
#include "mp.c"

static uint64_t mod[MP_SIZE] = { 0xffffffffffffffff, MWC_A1 - 1 };


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int m = MWC_A1 * b - 1;
   cpp_int r = cpp_int("0x2f65fed2e8400983a72f9a3547208003");
   cpp_int s = ((x + c * b) * r) % m;
   x = uint64_t(s);
   c = uint64_t(s >> 64);
*/

void jump(void) {
	static uint64_t jump[MP_SIZE] = { 0xa72f9a3547208003, 0x2f65fed2e8400983 };
	uint64_t state[MP_SIZE] = { x, c };
	mul(state, jump, mod);
	x = state[0];
	c = state[1];
}


/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations. 
   
   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int m = MWC_A1 * b - 1;
   cpp_int r = cpp_int("0x394649cfd6769c91e6f7814467f3fcdd");
   cpp_int s = ((x + c * b) * r) % m;
   x = uint64_t(s);
   c = uint64_t(s >> 64);
*/

void long_jump(void) {
	static uint64_t long_jump[MP_SIZE] = { 0xe6f7814467f3fcdd, 0x394649cfd6769c91 };
	uint64_t state[MP_SIZE] = { x, c };
	mul(state, long_jump, mod);
	x = state[0];
	c = state[1];
}