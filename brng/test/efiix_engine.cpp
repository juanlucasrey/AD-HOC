#include "../../test_simple/test_simple_include.hpp"
#include "seed_seq_inserter.hpp"
#include "test_tools_rng.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#include <efiix_engine.hpp>
#include <tools/mask.hpp>

#include <chrono>
#include <iostream>

#include "read_csv.hpp"

auto main() -> int {
    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint64_t>("./randomgen/efiix_state.txt"));
        adhoc::efiix64x48 rng(seq);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/efiix_vals.txt");

        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = values_from_python[i];
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::efiix64x48 rng1(1234);
        adhoc::efiix64x48 rng2(1234);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::efiix8x48 rng1(1);
        adhoc::efiix8x48 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw8();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::efiix16x48 rng1(1);
        adhoc::efiix16x48 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw16();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::efiix32x48 rng1(1);
        adhoc::efiix32x48 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::efiix64x48 rng1(1);
        adhoc::efiix64x48 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::efiix32x48 rng(1);
        check_fwd_and_back(rng, 1000000);
        adhoc::efiix32x48 rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::efiix32x48 rng(1);
        check_back_and_fwd(rng, 1000000);
        adhoc::efiix32x48 rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::efiix_engine<std::uint8_t, 8, 3> rng1(1);
        adhoc::efiix_engine<std::uint16_t, 8, 3> rng2(1);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::efiix_engine<std::uint32_t, 32, 13> rng1(1);
        adhoc::efiix_engine<std::uint64_t, 32, 13> rng2(1);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::efiix_engine<std::uint64_t, 32, 13> rng(1);
        check_fwd_and_back(rng, 1000000);
        adhoc::efiix_engine<std::uint64_t, 32, 13> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::efiix64x48>);

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::efiix64x48 rng(1234);
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "fwd" << std::endl;
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
            std::cout << "back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res1b);
    }

    TEST_END;
}
