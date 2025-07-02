#include "taus88.h"
#include <stdint.h>

#define SEED 987654321
static uint32_t s1 = SEED, s2 = SEED, s3 = SEED, b;
double taus88() {
    /* Generates numbers between 0 and 1. */
    b = (((s1 << 13) ^ s1) >> 19);
    s1 = (((s1 & 4294967294) << 12) ^ b);
    b = (((s2 << 2) ^ s2) >> 25);
    s2 = (((s2 & 4294967288) << 4) ^ b);
    b = (((s3 << 3) ^ s3) >> 11);
    s3 = (((s3 & 4294967280) << 17) ^ b);
    return ((s1 ^ s2 ^ s3) * 2.3283064365386963e-10);
}
