#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <jsf_engine.hpp>

#include "third_party/imneme/jsf.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#include "read_csv.hpp"

auto main() -> int {

    {
        auto const values_model =
            adhoc::readCSV<std::uint64_t>("./randomgen/jsf_state.txt");
        adhoc::jsf64 rng(values_model[0], values_model[1], values_model[2],
                         values_model[3], 0);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/jsf_vals.txt");
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        jsf32na rng1;
        adhoc::jsf32na rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        PractRand::RNGs::LightWeight::arbee rng1(1);
        adhoc::arbee rng2(1U, 1U, 2U, 3U, 12);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::jsf32 rng1(1);
        adhoc::jsf32 rng2(4058668781U, 1U, 1U, 1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::jsf64 rng1(1);
        adhoc::jsf64na rng2(4058668781U, 1U, 1U, 1U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::arbee rng(1U, 1U, 2U, 3U, 12);
        check_fwd_and_back(rng, 1000000);
        adhoc::arbee rng2(1U, 1U, 2U, 3U, 12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::arbee rng(1U, 1U, 2U, 3U, 12);
        check_back_and_fwd(rng, 1000000);
        adhoc::arbee rng2(1U, 1U, 2U, 3U, 12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::jsf32na rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::jsf32na rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::jsf32na rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::jsf32na rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::jsf32na rng1;
        adhoc::jsf_engine<std::uint64_t, 32, 27, 17, 0> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::jsf_engine<std::uint64_t, 32, 27, 17, 0> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::jsf_engine<std::uint64_t, 32, 27, 17, 0> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::jsf_engine<std::uint64_t, 32, 27, 17, 0> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::jsf_engine<std::uint64_t, 32, 27, 17, 0> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::jsf32>);

    TEST_END;
}
