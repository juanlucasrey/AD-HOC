#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#include <chrono>
#include <iostream>

#include <trivium_engine.hpp>

auto main() -> int {
    {
        PractRand::RNGs::LightWeight::trivium rng1(1);
        adhoc::trivium_engine<std::uint64_t> rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::trivium_engine<std::uint64_t> rng(1);
        check_fwd_and_back(rng, 1000000);
        adhoc::trivium_engine<std::uint64_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::trivium_engine<std::uint64_t> rng(1);
        check_back_and_fwd(rng, 1000000);
        adhoc::trivium_engine<std::uint64_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::trivium_engine<std::uint64_t> rng(1);
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

        std::uint64_t res2 = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res2 += rng.operator()<false>();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res2);
    }

    TEST_END;
}
