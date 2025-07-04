#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/imneme/xoroshiro.hpp"
#include <tools/uint128.hpp>
#include <xoroshiro_engine.hpp>

#include "read_csv.hpp"
#include "seed_seq_inserter.hpp"

auto main() -> int {

    {
        adhoc::seed_seq_inserter seq(adhoc::readCSV<std::uint64_t>(
            "./randomgen/xoroshiro128_state.txt"));
        adhoc::xoroshiro128plus rng(seq);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/xoroshiro128_vals.txt");
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::seed_seq_inserter seq(adhoc::readCSV<std::uint64_t>(
            "./randomgen/xoroshiro128pp_state.txt"));
        adhoc::xoroshiro128plusplus rng(seq);

        auto const values_from_python = adhoc::readCSV<std::uint64_t>(
            "./randomgen/xoroshiro128pp_vals.txt");
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        xoroshiro128plus64v0_1 rng1;
        adhoc::xoroshiro_engine<std::uint64_t, 64, 128, 55, 14, 36> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro128plus64v1_0 rng1;
        adhoc::xoroshiro128plus rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro128starstar64v1_0 rng1;
        adhoc::xoroshiro128starstar rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

#ifndef _MSC_VER
    {
        xoroshiro256plus128 rng1;
        adhoc::xoroshiro_engine<adhoc::uint128, 128, 256, 105, 36, 70> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }
#endif

    {
        xoroshiro64plus32vMEO1 rng1;
        adhoc::xoroshiro_engine<std::uint32_t, 32, 64, 27, 7, 20> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro64plus32v1_0 rng1;
        adhoc::xoroshiro_engine<std::uint32_t, 32, 64, 26, 9, 13> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro64star32v1_0 rng1;
        adhoc::xoroshiro64star rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro64starstar32v1_0 rng1;
        adhoc::xoroshiro64starstar rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro32plus16 rng1;
        adhoc::xoroshiro_engine<std::uint16_t, 16, 32, 13, 5, 10> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro32star16 rng1;
        adhoc::xoroshiro_engine<std::uint16_t, 16, 32, 13, 5, 10, 0x79BB> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro16plus8 rng1;
        adhoc::xoroshiro_engine<std::uint8_t, 8, 16, 4, 7, 3> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoroshiro16star8 rng1;
        adhoc::xoroshiro_engine<std::uint8_t, 8, 16, 4, 7, 3, 0xBB> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoroshiro_engine<std::uint32_t, 32, 64, 26, 9, 13> rng1;
        adhoc::xoroshiro_engine<std::uint64_t, 32, 64, 26, 9, 13> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoroshiro_engine<std::uint64_t, 32, 64, 26, 9, 13> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::xoroshiro_engine<std::uint64_t, 32, 64, 26, 9, 13> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoroshiro_engine<std::uint64_t, 32, 64, 26, 9, 13> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::xoroshiro_engine<std::uint64_t, 32, 64, 26, 9, 13> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoroshiro1024plusplus rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::xoroshiro1024plusplus rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoroshiro_engine<adhoc::uint128, 64, 128, 24, 16, 37, 5, 7, 9>
            rng1;
        adhoc::xoroshiro128starstar rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoroshiro_engine<std::uint64_t, 32, 64, 26, 9, 13, 0x9E3779BB>
            rng1;
        adhoc::xoroshiro64star rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoroshiro1024plusplus rng1;
        adhoc::xoroshiro_engine<adhoc::uint128, 64, 1024, 25, 27, 36, 0, 23, 0,
                                true>
            rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoroshiro128plusplus rng1;
        adhoc::xoroshiro_engine<adhoc::uint128, 64, 128, 49, 21, 28, 0, 17>
            rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::xoroshiro128plus>);

    TEST_END;
}
