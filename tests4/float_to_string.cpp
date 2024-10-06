///
/// Author: Daniel Wright
/// Influenced by:
/// https://github.com/vinceallenvince/floatToString/blob/master/floatToString.h
///
/// Simple demo involving a custom float to string implementation
/// Should be portable, but cannot process IEEE-754 double precision numbers
/// ^^ this can only handle single (float32) floating point numbers.
///
/// Decently accurate compared to sprintf, but I don't think it is going to be
/// nearly as rigorous.
///
/// Appears to be about ~5 times faster than glibc's sprintf("%f", number)
/// implementation, but that is at least on the machine tested. The performance
/// advantage seems to diminish slightly for numbers
/// with a very large magnitude (> (+/-)10^12)
///
/// compiled with following command:
/// gcc -O3 -march=native -malign-data=cache -fomit-frame-pointer
/// floatToString.cpp -o floatToString -lstdc++
///
/// used gcc version 6.1.1 on Arch Linux (Linux 4.7.0 (Intel x64))
///

#include <stdint.h>
#include <stdio.h>

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <gtest/gtest.h>

/// Random huge number to test how much divergence from
/// glibc's sprintf the implementation will acheieve
#define HI 1.6212959e+17
#define LO -HI

/// 250 thousand iterations, takes about a few seconds.
/// Covers a good amount of ground computationally for
/// just a simple benchmark.
#define ITERATIONS_PERF 250000

/// Used in function below to print side-by-side outputs
/// from both libc's sprintf and the custom implementation
#define ITERATIONS_LOGICAL 33

/// Used for char buffer
#define BUFFERSIZE 4096

/// Change to refine precision
#define NUM_DECIMAL_PLACES 15

/// ^
#define SPRINTF_STR "%15.15f"

/// Our ***very*** constrained digit-to-char function
static inline void to_digit_char_unsafe(int8_t digit, char *out) {
    out[0] = '0' + digit;
}

/// lossy floating point values probably not a big deal.
/// Enabled it since some platforms (e.g. ARM) can't use
/// some SIMD instructions (i.e. NEON) without relaxed fp.
/// Should only affect infinities, but I'm not entirely sure.
#pragma FP_CONTRACT fast
#define __FAST_MATH__ 1

/// Instead of calculating powers of 10 (and their reciprocals),
/// it is computationally quicker to have tables of these values
/// to which the code can look. Only goes up to 10^38 (and 10^-38)
/// since that is the order of magnitude of the maximum value for
/// an IEEE-754 32 bit floating point number.
/// See: https://en.wikipedia.org/wiki/Single-precision_floating-point_format
///
/// Note that a double (64-bit) has a maximum on the order of 10^308.
/// Feel free to implement that, but I didn't really have a demand for
/// implementing this for double precision as well
/// (float is enough precision for my use).
static const float __inv_pow_10[] = {
    1.0f,                                      // 0
    0.1f,                                      // 1
    0.01f,                                     // 2
    0.001f,                                    // 3
    0.0001f,                                   // 4
    0.00001f,                                  // 5
    0.000001f,                                 // 6
    0.0000001f,                                // 7
    0.00000001f,                               // 8
    0.000000001f,                              // 9
    0.0000000001f,                             // 10
    0.00000000001f,                            // 11
    0.000000000001f,                           // 12
    0.0000000000001f,                          // 13
    0.00000000000001f,                         // 14
    0.000000000000001f,                        // 15
    0.0000000000000001f,                       // 16
    0.00000000000000001f,                      // 17
    0.000000000000000001f,                     // 18
    0.0000000000000000001f,                    // 19
    0.00000000000000000001f,                   // 20
    0.000000000000000000001f,                  // 21
    0.0000000000000000000001f,                 // 22
    0.00000000000000000000001f,                // 23
    0.000000000000000000000001f,               // 24
    0.0000000000000000000000001f,              // 25
    0.00000000000000000000000001f,             // 26
    0.000000000000000000000000001f,            // 27
    0.0000000000000000000000000001f,           // 28
    0.00000000000000000000000000001f,          // 29
    0.000000000000000000000000000001f,         // 30
    0.0000000000000000000000000000001f,        // 31
    0.00000000000000000000000000000001f,       // 32
    0.000000000000000000000000000000001f,      // 33
    0.0000000000000000000000000000000001f,     // 34
    0.00000000000000000000000000000000001f,    // 35
    0.000000000000000000000000000000000001f,   // 36
    0.00000000000000000000000000000000000001f, // 37
    0.000000000000000000000000000000000000001f // 38

};

/// Same concept as above, but positive exponents of 10.
static const float __pos_pow_10[] = {

    1.0f,                                      // 0
    10.0f,                                     // 1
    100.0f,                                    // 2
    1000.0f,                                   // 3
    10000.0f,                                  // 4
    100000.0f,                                 // 5
    1000000.0f,                                // 6
    10000000.0f,                               // 7
    100000000.0f,                              // 8
    1000000000.0f,                             // 9
    10000000000.0f,                            // 10
    100000000000.0f,                           // 11
    1000000000000.0f,                          // 12
    10000000000000.0f,                         // 13
    100000000000000.0f,                        // 14
    1000000000000000.0f,                       // 15
    10000000000000000.0f,                      // 16
    100000000000000000.0f,                     // 17
    1000000000000000000.0f,                    // 18
    10000000000000000000.0f,                   // 19
    100000000000000000000.0f,                  // 20
    1000000000000000000000.0f,                 // 21
    10000000000000000000000.0f,                // 22
    100000000000000000000000.0f,               // 23
    1000000000000000000000000.0f,              // 24
    10000000000000000000000000.0f,             // 25
    100000000000000000000000000.0f,            // 26
    1000000000000000000000000000.0f,           // 27
    10000000000000000000000000000.0f,          // 28
    100000000000000000000000000000.0f,         // 29
    1000000000000000000000000000000.0f,        // 30
    10000000000000000000000000000000.0f,       // 31
    100000000000000000000000000000000.0f,      // 32
    1000000000000000000000000000000000.0f,     // 33
    10000000000000000000000000000000000.0f,    // 34
    100000000000000000000000000000000000.0f,   // 35
    1000000000000000000000000000000000000.0f,  // 36
    10000000000000000000000000000000000000.0f, // 37
    100000000000000000000000000000000000000.0f // 38

};

/// Offset pointers for calculation below (to save a few instructions)
static const float *const __inv_pow_10p = &(__inv_pow_10[-1]);
static const float *const __pos_pow_10o = &(__pos_pow_10[1]);
static const float *const __pos_pow_10p = &(__pos_pow_10[-1]);

static size_t floatToString(char *outstr, float value,
                            size_t places = NUM_DECIMAL_PLACES) {

    /// Used to store current digit for calculations below.
    /// That means the possible values for this byte are [0, 9]
    int8_t digit;

    /// The order of magnitude on which the input resides
    /// i.e. tenscount = (int)log10f(abs(value))
    /// ^^ dont use log10f here -- tried and it's miserably slow :-)
    size_t tenscount = 0;

    /// counter for loops below
    size_t i;

    /// temp variable for calculations below
    float tempfloat = value;

    /// the number of characters consumed by the string representing this float
    /// (this is what is returned)
    size_t c = 0;

    /// calculate rounding term d:   0.5/pow(10,places)
    float d = (value < 0) ? -0.5f : 0.5f;

    /// divide by ten for each decimal place
    d = d * __inv_pow_10[places];

    /// this small addition, combined with truncation will round our values
    /// properly
    if (value < 0) {
        tempfloat = -(value + d);
    } else {
        tempfloat = value + d;
    }

    /// This logic is similar to
    /// tenscount = (int)log10f(tempfloat)
    /// but is not miserably slow.
    while (__pos_pow_10o[tenscount++] <= tempfloat) {
    }

    /// the number is negative
    if (value < 0) {
        outstr[c++] = '-';
    }

    if (tenscount == 0) {
        outstr[c++] = '0';
    }

#pragma nounroll
    for (i = 0; i < tenscount; i++) {

        const ssize_t idx = tenscount - i;

        digit = (int8_t)(tempfloat * __inv_pow_10p[idx]);

        to_digit_char_unsafe(digit, &outstr[c++]);
        tempfloat = tempfloat - ((float)digit * __pos_pow_10p[idx]);
    }

    /// if places is zero, then there is no decimal part
    if (places > 0) {
        outstr[c++] = '.';
    }

#pragma nounroll
    for (i = 0; i < places; i++) {

        tempfloat *= 10.0f;
        digit = (int8_t)tempfloat;

        /// convert digit to character
        /// If you replace this function with sprintf("%d", digit)
        /// you will see where sprintf most likely spends most of
        /// its time during their implementation of this routine :-)
        ///
        /// Note that if for some reason digit is not within [0,9],
        /// your number will be a flaming pile of shit
        to_digit_char_unsafe(digit, &outstr[c++]);

        tempfloat = tempfloat - (float)digit;
    }

    outstr[c++] = '\0';

    return c;
}

static inline float getRandom() {
    return static_cast<float>(LO) +
           static_cast<float>(rand()) /
               (static_cast<float>(RAND_MAX / (HI - LO)));
}

static int perf_test() {

    int SZ = ITERATIONS_PERF;

    clock_t t1, t2;

    double nTime, cTime, startTime, totalTime;

    nTime = cTime = 0.0;

    char bfr[BUFFERSIZE];

    const size_t max_n = BUFFERSIZE - 1;

    /// set seed for our random
    srand(static_cast<unsigned>(time(0)));

    float *randoms = new float[SZ];

    for (int i = 0; i < SZ; i++) {
        randoms[i] = getRandom();
    }

    ////////////////////////////////
    /// sprintf

    t1 = clock();

    for (int i = 0; i < SZ; i++) {
        snprintf(bfr, max_n, SPRINTF_STR, static_cast<double>(randoms[i]));
    }

    t2 = clock();

    startTime = (double)t1;
    nTime = ((double)(t2 - t1) / 1000000.0) * 1000.0;

    ////////////////////////////////
    /// custom format

    t1 = clock();

    for (int i = 0; i < SZ; i++) {
        floatToString(bfr, randoms[i]);
    }

    t2 = clock();

    cTime = ((double)(t2 - t1) / 1000000.0) * 1000.0;
    // totalTime = ((double)(t2 - startTime) / 1000000.0) * 1000.0;

    // printf("processed %d random numbers\n", SZ);

    // printf("Done processing."
    //        "\nTotal time spent: %6.6f ms"
    //        "\nsprintf		: %6.6f ms"
    //        "\ncustom		: %6.6f ms\n",
    //        totalTime, nTime, cTime);

    if (nTime < cTime) {
        printf("native sprintf is %6.6f times faster\n", (cTime / nTime));
    } else if (nTime > cTime) {
        printf("custom is %6.6f times faster\n", (nTime / cTime));
    } else {
        printf("Both exhibited the same exact performance.\n");
    }

    delete[] randoms;

    return 0;
}

static int logical_test() {

    char bfr[BUFFERSIZE];
    const size_t max_n = BUFFERSIZE - 1;

    /// set seed for "random"
    srand(static_cast<unsigned>(time(0)));

    int SZ = ITERATIONS_LOGICAL;

    float *randoms = new float[SZ];

    for (int i = 0; i < SZ; i++) {
        randoms[i] = getRandom();
    }

    printf("sprintf impl\t\tfloatToString impl\n\n");

    for (int i = 0; i < SZ; i++) {

        snprintf(bfr, max_n, SPRINTF_STR, static_cast<double>(randoms[i]));
        printf("%s\t", bfr);

        floatToString(bfr, randoms[i]);
        printf("%s\n", bfr);
    }

    delete[] randoms;

    return 0;
}

// int main() {
TEST(Sprintf, Init2) {

    // perf_test();

    // printf("\n");

    // logical_test();
    float test = 9.9905001000000044e+30f;
    char bfr[BUFFERSIZE];
    floatToString(bfr, test);
    std::cout << bfr << std::endl;

    // return 0;
}
