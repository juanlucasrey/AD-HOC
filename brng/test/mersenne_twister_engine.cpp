#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "../include/mersenne_twister_engine.hpp"
#include "seed_seq_inserter.hpp"

#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#define REAL_GCC __GNUC__ // probably
#endif

auto main() -> int {

    {
        std::seed_seq seq1{1, 2, 3, 4, 5};
        adhoc::mt19937_64 rng1(seq1);
        std::seed_seq seq2{1, 2, 3, 4, 5};
        std::mt19937_64 rng2(seq2);
        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    // std check
    {
        adhoc::mt19937 gen32;    // overload (1)
        adhoc::mt19937_64 gen64; // overload (1)
        gen32.discard(10000 - 1);
        gen64.discard(10000 - 1);
        EXPECT_EQUAL(gen32(), 4123659995);
        EXPECT_EQUAL(gen64(), 9981545732273789042ull);
    }

    {
        static_assert(adhoc::mt19937::min() == 0U);
        static_assert(adhoc::mt19937::max() == 4294967295U);
    }

    // check against std
    {
        std::mt19937 rng1;
        adhoc::mt19937 rng2;
        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        PractRand::RNGs::LightWeight::mt19937 rng1(5489U);
        adhoc::mt19937 rng2;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    // check against std 64
    {
        std::mt19937_64 rng1;
        adhoc::mt19937_64 rng2;
        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    // std flaw 1
    {
        {
            std::vector<std::uint32_t> seed;
            seed.resize(624);
            std::iota(seed.begin(), seed.end(), 2);
            adhoc::seed_seq_inserter seq(seed);

            std::mt19937 rng1(seq);

            seed.front()++;
            adhoc::seed_seq_inserter seq2(seed);
            std::mt19937 rng2(seq2);
            // both rng are NOT seen as equal initially
            EXPECT_NOT_EQUAL(rng1, rng2);

            // however the values they output are equal!:
            // 1- first value is equal
            EXPECT_EQUAL(rng1(), rng2());
            // 2- after first value is output, rng are now seen as equal
            EXPECT_EQUAL(rng1, rng2);
            // and of course they generate same values
            compare_rng(rng1, rng2, 624);
        }

        // this does not happen with adhoc
        {
            std::vector<std::uint32_t> seed;
            seed.resize(624);
            std::iota(seed.begin(), seed.end(), 2);
            adhoc::seed_seq_inserter seq(seed);

            adhoc::mt19937 rng1(seq);

            seed.front()++;
            adhoc::seed_seq_inserter seq2(seed);

            adhoc::mt19937 rng2(seq2);
            // both rng are equal
            EXPECT_EQUAL(rng1, rng2);
            // and of course they generate same values
            compare_rng(rng1, rng2, 624);
        }
    }

    // std flaw 2 (only in gcc)
    {
#if defined(REAL_GCC)
        {
            std::vector<std::uint32_t> seed;
            seed.resize(624);
            std::iota(seed.begin(), seed.end(), 2);
            adhoc::seed_seq_inserter seq(seed);
            std::mt19937 rng1(seq);
            rng1.discard(1);

            std::iota(seq.vals.begin(), seq.vals.end() - 1, 3);
            seq.vals.back() = 2567483729U;
            std::mt19937 rng2(seq);

            // both rng are not seen as equal
            EXPECT_NOT_EQUAL(rng1, rng2);

            // but they generate same values! (624 values is enough to
            // deteremine they are equal)
            compare_rng(rng1, rng2, 624);

            // same check again..
            EXPECT_NOT_EQUAL(rng1, rng2);
            compare_rng(rng1, rng2, 624);
        }
#endif

// this does not happen with windows cl
// and SOMETIMES with CLANG
#ifdef _MSC_VER
        // this does not happen with clang and windows cl
        {
            std::vector<std::uint32_t> seed;
            seed.resize(624);
            std::iota(seed.begin(), seed.end(), 2);
            adhoc::seed_seq_inserter seq(seed);
            std::mt19937 rng1(seq);
            auto val = rng1();

            std::iota(seq.vals.begin(), seq.vals.end() - 1, 3);
            seq.vals.back() = 2567483729U;
            std::mt19937 rng2(seq);

            // both rng are seen as equal
            EXPECT_EQUAL(rng1, rng2);

            // and they of course generate the same values
            compare_rng(rng1, rng2, 624);
        }
#endif

        // this does not happen with adhoc
        {
            std::vector<std::uint32_t> seed;
            seed.resize(624);
            std::iota(seed.begin(), seed.end(), 2);
            adhoc::seed_seq_inserter seq(seed);

            adhoc::mt19937 rng1(seq);
            rng1();

            std::iota(seed.begin(), seed.end() - 1, 3);
            seed.back() = 2567483729U;
            adhoc::seed_seq_inserter seq2(seed);
            adhoc::mt19937 rng2(seq2);

            // both rng are seen as equal
            EXPECT_EQUAL(rng1, rng2);

            // and they of course generate the same values
            compare_rng(rng1, rng2, 624);
        }
    }

    {
        adhoc::mt19937 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::mt19937 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mt19937 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::mt19937 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mt19937_64 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::mt19937_64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }
    {
        adhoc::mt19937_64 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::mt19937_64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    // 32 many types
    {
        adhoc::mersenne_twister_engine<
            std::uint32_t, 32, 624, 397, 31, 0x9908b0df, 11, 0xffffffff, 7,
            0x9d2c5680, 15, 0xefc60000, 18, 1812433253>
            rng1(1234), rng1b(1234);
        adhoc::mersenne_twister_engine<
            std::uint64_t, 32, 624, 397, 31, 0x9908b0df, 11, 0xffffffff, 7,
            0x9d2c5680, 15, 0xefc60000, 18, 1812433253>
            rng2(1234), rng2b(1234);

        check_back_and_fwd(rng1, 1000000);
        check_back_and_fwd(rng2, 1000000);

        EXPECT_EQUAL(rng1, rng1b);
        EXPECT_EQUAL(rng2, rng2b);

#ifndef _MSC_VER
        adhoc::mersenne_twister_engine<
            __uint128_t, 32, 624, 397, 31, 0x9908b0df, 11, 0xffffffff, 7,
            0x9d2c5680, 15, 0xefc60000, 18, 1812433253>
            rng3(1234), rng3b(1234);
        check_back_and_fwd(rng3, 1000000);
        EXPECT_EQUAL(rng3, rng3b);
#endif

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);

#ifndef _MSC_VER
            auto val3 = rng3b();
            EXPECT_EQUAL(val1, val3);
#endif
        }
    }

#ifndef _MSC_VER
    // 64 many types
    {
        adhoc::mersenne_twister_engine<
            std::uint64_t, 64, 312, 156, 31, 0xb5026f5aa96619e9, 29,
            0x5555555555555555, 17, 0x71d67fffeda60000, 37, 0xfff7eee000000000,
            43, 6364136223846793005>
            rng1;
        adhoc::mersenne_twister_engine<
            __uint128_t, 64, 312, 156, 31, 0xb5026f5aa96619e9, 29,
            0x5555555555555555, 17, 0x71d67fffeda60000, 37, 0xfff7eee000000000,
            43, 6364136223846793005>
            rng2;

        check_back_and_fwd(rng1, 1000000);
        check_back_and_fwd(rng2, 1000000);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }
#endif

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::mt19937 rng;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new" << std::endl;
            std::cout << time << std::endl;
        }

        std::uint64_t res2 = 0;
        std::mt19937 rng2;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res2 += rng2();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "old" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res2);

        std::uint64_t res1b = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1b += rng.operator()<false>();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res1b);
    }

    TEST_END;
}
