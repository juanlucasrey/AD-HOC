#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/imneme/xoshiro.hpp"
#include <xoshiro_engine.hpp>

#include "read_csv.hpp"
#include "seed_seq_inserter.hpp"

auto main() -> int {

    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint64_t>("./randomgen/xoshiro256_state.txt"));
        adhoc::xoshiro256starstar64 rng(seq);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/xoshiro256_vals.txt");
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint64_t>("./randomgen/xoshiro512_state.txt"));
        adhoc::xoshiro512starstar64 rng(seq);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/xoshiro512_vals.txt");
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        xoshiro512plus64v1_0 rng1;
        adhoc::xoshiro512plus64v1_0 rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro512starstar64v1_0 rng1;
        adhoc::xoshiro512starstar64v1_0 rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro256plus64v1_0 rng1;
        adhoc::xoshiro256plus64v1_0 rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro256starstar64v1_0 rng1;
        adhoc::xoshiro256starstar64v1_0 rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro128plus32v1_0 rng1;
        adhoc::xoshiro128plus32v1_0 rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro128starstar32v1_0 rng1;
        adhoc::xoshiro128starstar32v1_0 rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro64plus16xxx rng1;
        adhoc::xoshiro64plus16xxx rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro64starstar16xxx rng1;
        adhoc::xoshiro64starstar16xxx rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro32plus8xxx rng1;
        adhoc::xoshiro32plus8xxx rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro32plus8yyy rng1;
        adhoc::xoshiro32plus8yyy rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro32starstar8xxx rng1;
        adhoc::xoshiro32starstar8xxx rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        xoshiro32starstar8yyy rng1;
        adhoc::xoshiro32starstar8yyy rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoshiro128plus32 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::xoshiro128plus32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro128plus32 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::xoshiro128plus32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro512plus64 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::xoshiro512plus64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro512plus64 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::xoshiro512plus64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint32_t, 32, 4, 9, 11> rng1;
        adhoc::xoshiro_engine<std::uint64_t, 32, 4, 9, 11> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint64_t, 32, 4, 9, 11> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::xoshiro_engine<std::uint64_t, 32, 4, 9, 11> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint64_t, 32, 4, 9, 11> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::xoshiro_engine<std::uint64_t, 32, 4, 9, 11> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint32_t, 32, 4, 9, 11, 5, 7, 9> rng1;
        adhoc::xoshiro_engine<std::uint64_t, 32, 4, 9, 11, 5, 7, 9> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint32_t, 32, 4, 9, 11, 5, 7, 9> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::xoshiro_engine<std::uint32_t, 32, 4, 9, 11, 5, 7, 9> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint32_t, 32, 4, 9, 11, 5, 7, 9> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::xoshiro_engine<std::uint32_t, 32, 4, 9, 11, 5, 7, 9> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint32_t, 32, 8, 9, 11, 5, 7, 9> rng1;
        adhoc::xoshiro_engine<std::uint64_t, 32, 8, 9, 11, 5, 7, 9> rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint32_t, 32, 8, 9, 11, 5, 7, 9> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::xoshiro_engine<std::uint32_t, 32, 8, 9, 11, 5, 7, 9> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::xoshiro_engine<std::uint32_t, 32, 8, 9, 11, 5, 7, 9> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::xoshiro_engine<std::uint32_t, 32, 8, 9, 11, 5, 7, 9> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::xoshiro256starstar64>);
    static_assert(std::bidirectional_iterator<adhoc::xoshiro512starstar64>);

    TEST_END;
}
