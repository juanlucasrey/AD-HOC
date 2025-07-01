#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#include <chrono>
#include <iostream>

#include <isaac_engine.hpp>
#include <tools/uint128.hpp>

extern "C" {
#include "external/burtleburtle/isaac.h"
}

auto main() -> int {

    {
        PractRand::RNGs::LightWeight::isaac32x256 rng1(1);
        adhoc::isaac32x256 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::isaac64x256 rng1(1);
        adhoc::isaac64x256 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::isaac32x256 rng(1);
        check_fwd_and_back(rng, 10000000);
        adhoc::isaac32x256 rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::isaac32x256 rng(1);
        check_back_and_fwd(rng, 10000000);
        adhoc::isaac32x256 rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::isaac64x256 rng(1);
        check_fwd_and_back(rng, 10000000);
        adhoc::isaac64x256 rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::isaac64x256 rng(1);
        check_back_and_fwd(rng, 10000000);
        adhoc::isaac64x256 rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::isaac32x256>);
    static_assert(std::bidirectional_iterator<adhoc::isaac64x256>);

    {
        adhoc::isaac_engine<std::uint64_t, 64, 21, -5, 12, -33, 3, true,
                            0x9e3779b97f4a7c13ULL>
            rng1(12);

        adhoc::isaac_engine<adhoc::uint128, 64, 21, -5, 12, -33, 3, true,
                            0x9e3779b97f4a7c13ULL>
            rng2(12);

        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::isaac_engine<adhoc::uint128, 64, 21, -5, 12, -33, 3, true,
                            0x9e3779b97f4a7c13ULL>
            rng(12);
        check_fwd_and_back(rng, 1000000);
        adhoc::isaac_engine<adhoc::uint128, 64, 21, -5, 12, -33, 3, true,
                            0x9e3779b97f4a7c13ULL>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::isaac_engine<adhoc::uint128, 64, 21, -5, 12, -33, 3, true,
                            0x9e3779b97f4a7c13ULL>
            rng(12);
        check_back_and_fwd(rng, 1000000);
        adhoc::isaac_engine<adhoc::uint128, 64, 21, -5, 12, -33, 3, true,
                            0x9e3779b97f4a7c13ULL>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::isaac_engine<std::uint32_t, 32, 13, -6, 2, -16, 2, false,
                            0x9e3779b9UL>
            rng1(12);

        adhoc::isaac_engine<std::uint64_t, 32, 13, -6, 2, -16, 2, false,
                            0x9e3779b9UL>
            rng2(12);

        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::isaac_engine<std::uint64_t, 32, 13, -6, 2, -16, 2, false,
                            0x9e3779b9UL>
            rng(12);
        check_fwd_and_back(rng, 1000000);
        adhoc::isaac_engine<std::uint64_t, 32, 13, -6, 2, -16, 2, false,
                            0x9e3779b9UL>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::isaac_engine<std::uint64_t, 32, 13, -6, 2, -16, 2, false,
                            0x9e3779b9UL>
            rng(12);
        check_back_and_fwd(rng, 1000000);
        adhoc::isaac_engine<std::uint64_t, 32, 13, -6, 2, -16, 2, false,
                            0x9e3779b9UL>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint32_t res1 = 0;
        adhoc::isaac32x256 rng(std::uint32_t(1234));
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                auto const val = rng();
                res1 += val > 0 ? val : (val + 1.0);
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "fwd" << std::endl;
            std::cout << time << std::endl;
        }

        std::uint32_t res1b = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                auto const val = rng.operator()<false>();
                res1b += val > 0 ? val : (val + 1.0);
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res1b);
    }

    TEST_END;
}
