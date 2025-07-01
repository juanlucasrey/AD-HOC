#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <tools/uint128.hpp>
#include <xsm_engine.hpp>

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#include <chrono>

auto main() -> int {

    {
        PractRand::RNGs::LightWeight::xsm64 rng1(1);
        adhoc::xsm64 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::xsm32 rng1(1);
        adhoc::xsm32 rng2(1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::xsm32 rng(1U);
        check_fwd_and_back(rng, 2);
        adhoc::xsm32 rng2(1U);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xsm32 rng(1U);
        check_back_and_fwd(rng, 1000000);
        adhoc::xsm32 rng2(1U);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xsm64 rng(1U);
        check_fwd_and_back(rng, 1000000);
        adhoc::xsm64 rng2(1U);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xsm64 rng(1U);
        check_back_and_fwd(rng, 1000000);
        adhoc::xsm64 rng2(1U);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xsm_engine<std::uint32_t, 32, 9, 19, 16, 0xD251CF2D, 0x299529B5>
            rng1;
        adhoc::xsm_engine<std::uint64_t, 32, 9, 19, 16, 0xD251CF2D, 0x299529B5>
            rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xsm_engine<std::uint64_t, 32, 9, 19, 16, 0xD251CF2D, 0x299529B5>
            rng(1U);
        check_fwd_and_back(rng, 1000000);
        adhoc::xsm_engine<std::uint64_t, 32, 9, 19, 16, 0xD251CF2D, 0x299529B5>
            rng2(1U);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xsm_engine<std::uint_fast64_t, 64, 16, 40, 32,
                          0xA3EC647659359ACDULL>
            rng1;
        adhoc::xsm_engine<adhoc::uint128, 64, 16, 40, 32, 0xA3EC647659359ACDULL>
            rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xsm_engine<adhoc::uint128, 64, 16, 40, 32, 0xA3EC647659359ACDULL>
            rng(1U);
        check_fwd_and_back(rng, 1000000);
        adhoc::xsm_engine<adhoc::uint128, 64, 16, 40, 32, 0xA3EC647659359ACDULL>
            rng2(1U);
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::xsm32>);
    static_assert(std::bidirectional_iterator<adhoc::xsm64>);

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint32_t res1 = 0;
        adhoc::xsm64 rng(std::uint32_t(1234));
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
