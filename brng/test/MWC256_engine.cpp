#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/unimi/MWC256.c"
#include <mwc_engine.hpp>
#include <tools/mask.hpp>
#include <tools/modular_arithmetic.hpp>

auto main() -> int {

#ifndef _MSC_VER
    {
        adhoc::mcw256 rng(9012, 5678, 1234, 1);

        z = 9012;
        y = 5678;
        x = 1234;
        c = 1;
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = next();
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }
#endif

    {
        adhoc::mcw256 rng(1234, 5678, 9012, 1);
        check_fwd_and_back(rng, 1000000);
        adhoc::mcw256 rng2(1234, 5678, 9012, 1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcw256 rng(1234, 5678, 9012, 1);
        check_back_and_fwd(rng, 1000000);
        adhoc::mcw256 rng2(1234, 5678, 9012, 1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcw_engine<std::uint32_t, 32, 128, 0xfff62cf2> rng1(2, 3, 4, 5);
        adhoc::mcw_engine<std::uint64_t, 32, 128, 0xfff62cf2> rng2(2, 3, 4, 5);
        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::mcw_engine<std::uint64_t, 32, 128, 0xfff62cf2> rng(2, 3, 4, 5);
        check_fwd_and_back(rng, 1000000);
        adhoc::mcw_engine<std::uint64_t, 32, 128, 0xfff62cf2> rng2(2, 3, 4, 5);
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::mcw256>);

    TEST_END;
}
