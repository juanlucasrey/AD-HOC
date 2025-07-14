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

/*

This is a minimal, ludicrously inefficient library to perform modular
arithmetic on multiple-precision natural numbers represented by an
array of MP_SIZE uint64_t integers using 128-bit arithmetic.

The only purpose of this code is to make the jump code of (G)MWC
generators self contained. Using any proper multiprecision library
will be more efficient, albeit at this size the increase in speed
will not be so dramatic (in my tests, less than fourfold).

This code is intentially not documented to discourage general usage.

*/

#include <string.h>
#include <stdint.h>

static int cmp(const uint64_t * const a, const uint64_t * const b) {
	int i;
	for(i = MP_SIZE; i-- != 0; ) {
		if (a[i] < b[i]) return -1;
		if (a[i] > b[i]) return 1;
	}
	return 0;
}

// Assumes a >= b
static void _sub(uint64_t * const a, const uint64_t * const b) {
	int borrow = 0;
	for(int i = 0; i < MP_SIZE; i++) {
		__int128_t d = (__int128_t)a[i] - (__int128_t)b[i] - (__int128_t)borrow;
		borrow = d < 0;
		a[i] = ((__int128_t)UINT64_MAX + 1) + d;
	}
}

static void rem(uint64_t * const a, const uint64_t * const m) {
	for(;;) {
		if (cmp(a, m) < 0) return;
		_sub(a, m);
	}
}

static void add(uint64_t * const a, const uint64_t * const b, const uint64_t * const m) {
	int carry = 0;
	for(int i = 0; i < MP_SIZE; i++) {
		__uint128_t s = (__uint128_t)a[i] + (__uint128_t)b[i] + (__uint128_t)carry;
		carry = s > UINT64_MAX;
		a[i] = s;
	}
	if (m) rem(a, m);
}

static void sub(uint64_t * const a, const uint64_t * const b, const uint64_t * const m) {
	if (cmp(a, b) >= 0) _sub(a, b);
	else {
		uint64_t t[MP_SIZE];
		memcpy(t, m, sizeof t);
		_sub(t, b);
		add(a, t, m);
	}
}

static void mul(uint64_t * const a, const uint64_t * const b, const uint64_t * const m) {
	uint64_t r[MP_SIZE] = {}, t[MP_SIZE];
	memcpy(t, a, sizeof t);

	int d;
	for(d = MP_SIZE; d-- != 0 && b[d] == 0;);
	d++;
	for(int i = 0; i < d * 64; i++) {
		if (b[i >> 6] & (UINT64_C(1) << (i & 63))) add(r, t, m);
		add(t, t, m);
	}

	memcpy(a, r, sizeof r);
}