#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "../include/distribution/uniform_distribution.hpp"
#include "../include/lfsr_xor_engine.hpp"

#include "external/lfsr113.c"
#include "external/lfsr258.c"
#include "external/taus88.c"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>

auto main() -> int {

    {
        adhoc::canonical<adhoc::lfsr113<std::uint_fast32_t>> rng;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = lfsr113();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::canonical<adhoc::lfsr258<std::uint_fast64_t>> rng;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = lfsr258();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::canonical<adhoc::taus88<std::uint_fast64_t>> rng;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = taus88();
            EXPECT_EQUAL(val1, val2);
        }
    }

    // lfsr113
    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::lfsr113<std::uint32_t> rng(seq);
        check_fwd_and_back(rng, 1000000);
        adhoc::lfsr113<std::uint32_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::lfsr113<std::uint32_t> rng(seq);
        check_back_and_fwd(rng, 1000000);
        adhoc::lfsr113<std::uint32_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::lfsr113<std::uint64_t> rng(seq);
        check_fwd_and_back(rng, 1000000);
        adhoc::lfsr113<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::lfsr113<std::uint64_t> rng(seq);
        check_back_and_fwd(rng, 1000000);
        adhoc::lfsr113<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::lfsr113<std::uint32_t> rng1;
        adhoc::lfsr113<std::uint64_t> rng2;
        compare_rng(rng1, rng2, 1000000);
    }

    // lfsr258
    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::lfsr258<std::uint64_t> rng(seq);
        check_fwd_and_back(rng, 1000000);
        adhoc::lfsr258<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::lfsr258<std::uint64_t> rng(seq);
        check_back_and_fwd(rng, 1000000);
        adhoc::lfsr258<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

#ifndef _MSC_VER
    {
        adhoc::lfsr258<std::uint64_t> rng1;
        adhoc::lfsr258<__uint128_t> rng2;
        compare_rng(rng1, rng2, 1000000);
    }
#endif

    // taus88
    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::taus88<std::uint32_t> rng(seq);
        check_fwd_and_back(rng, 1000000);
        adhoc::taus88<std::uint32_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::taus88<std::uint32_t> rng(seq);
        check_back_and_fwd(rng, 1000000);
        adhoc::taus88<std::uint32_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::taus88<std::uint64_t> rng(seq);
        check_fwd_and_back(rng, 1000000);
        adhoc::taus88<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::taus88<std::uint64_t> rng(seq);
        check_back_and_fwd(rng, 1000000);
        adhoc::taus88<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::taus88<std::uint32_t> rng1;
        adhoc::taus88<std::uint64_t> rng2;
        compare_rng(rng1, rng2, 1000000);
    }

    int sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::lfsr113<std::uint_fast32_t> rng;
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

        adhoc::lfsr113<std::uint_fast64_t> rng64;
        std::uint64_t res64 = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res64 += rng64();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new 64" << std::endl;
            std::cout << time << std::endl;
        }

        EXPECT_EQUAL(res1, res64);
    }

    TEST_END;
}
