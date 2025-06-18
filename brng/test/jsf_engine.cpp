#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <jsf_engine.hpp>

#include "external/jsf.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

auto main() -> int {

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

    TEST_END;
}
