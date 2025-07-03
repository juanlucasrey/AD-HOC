#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <tyche_engine.hpp>

#include "read_csv.hpp"

#include "seed_seq_inserter.hpp"

auto main() -> int {
    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint32_t>("./randomgen/tyche_state.txt"));

        auto const values_from_python = adhoc::split_uint64_to_uint32(
            adhoc::readCSV<std::uint64_t>("./randomgen/tyche_vals.txt"), true);

        adhoc::tyche rng(seq);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::seed_seq_inserter seq(adhoc::readCSV<std::uint32_t>(
            "./randomgen/tyche_openrand_state.txt"));

        auto const values_from_python = adhoc::split_uint64_to_uint32(
            adhoc::readCSV<std::uint64_t>(
                "./randomgen/tyche_openrand_vals.txt"),
            true);

        adhoc::tyche_openrand rng(seq);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::tyche rng(123);
        check_fwd_and_back(rng, 1000000);
        adhoc::tyche rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::tyche rng(123);
        check_back_and_fwd(rng, 1000000);
        adhoc::tyche rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::tyche_openrand rng(123);
        check_fwd_and_back(rng, 1000000);
        adhoc::tyche_openrand rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::tyche_openrand rng(123);
        check_back_and_fwd(rng, 1000000);
        adhoc::tyche_openrand rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::tyche>);
    static_assert(std::bidirectional_iterator<adhoc::tyche_openrand>);

    {
        adhoc::tyche_engine<std::uint32_t, 32> rng1(12);
        adhoc::tyche_engine<std::uint64_t, 32> rng2(12);

        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::tyche_engine<std::uint32_t, 32, true> rng1(12);
        adhoc::tyche_engine<std::uint64_t, 32, true> rng2(12);

        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::tyche_engine<std::uint64_t, 32> rng(123);
        check_fwd_and_back(rng, 1000000);
        adhoc::tyche_engine<std::uint64_t, 32> rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::tyche_engine<std::uint64_t, 32> rng(123);
        check_back_and_fwd(rng, 1000000);
        adhoc::tyche_engine<std::uint64_t, 32> rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::tyche_engine<std::uint64_t, 32, true> rng(123);
        check_fwd_and_back(rng, 1000000);
        adhoc::tyche_engine<std::uint64_t, 32, true> rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::tyche_engine<std::uint64_t, 32, true> rng(123);
        check_back_and_fwd(rng, 1000000);
        adhoc::tyche_engine<std::uint64_t, 32, true> rng2(123);
        EXPECT_EQUAL(rng, rng2);
    }

    TEST_END;
}
