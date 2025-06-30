#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/imneme/sfc.hpp"
#include <sfc_engine.hpp>

#include "read_csv.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#include <chrono>
#include <iostream>

auto main() -> int {

    {
        auto const values_model =
            adhoc::readCSV<std::uint64_t>("./randomgen/sfc64_state.txt");
        adhoc::sfc64 rng(values_model[0], values_model[1], values_model[2],
                         values_model[3], 0);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/sfc64_vals.txt");
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::sfc16 rng1(1);
        adhoc::sfc16 rng2(1, 0, 0, 0, 10);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw16();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::sfc32 rng1(1);
        adhoc::sfc32 rng2(0, 1);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::sfc64 rng1(1);
        adhoc::sfc64 rng2(1, 1, 1);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        sfc64a rng1;
        adhoc::sfc64a rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc64b rng1;
        adhoc::sfc64b rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc32a rng1;
        adhoc::sfc32a rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc32b rng1;
        adhoc::sfc32b rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc32c rng1;
        adhoc::sfc32c rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc16a rng1;
        adhoc::sfc16a rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc16b rng1;
        adhoc::sfc16b rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc16c rng1;
        adhoc::sfc16c rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc16d rng1;
        adhoc::sfc16d rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc16e rng1;
        adhoc::sfc16e rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc16f rng1;
        adhoc::sfc16f rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        sfc8 rng1;
        adhoc::sfc8 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::sfc64a rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::sfc64a rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfc64a rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::sfc64a rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfc_engine<std::uint64_t, 32, 21, 9, 3> rng1;
        adhoc::sfc32 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::sfc_engine<std::uint64_t, 32, 21, 9, 3> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::sfc_engine<std::uint64_t, 32, 21, 9, 3> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::sfc64>);

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::sfc64 rng;
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

        std::uint64_t res1b = 0;
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
