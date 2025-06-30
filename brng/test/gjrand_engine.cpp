#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/imneme/gjrand.hpp"
#include <gjrand_engine.hpp>

auto main() -> int {

    {
        adhoc::gjrand32 rng1;
        adhoc::gjrand_engine<std::uint64_t, 32U, 0x96a5U, 16U, 11U, 19U> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        gjrand64 rng1;
        adhoc::gjrand64 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        gjrand32 rng1;
        adhoc::gjrand32 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        gjrand16 rng1;
        adhoc::gjrand16 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        gjrand8 rng1;
        adhoc::gjrand8 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::gjrand64 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::gjrand64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::gjrand64 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::gjrand64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::gjrand_engine<std::uint64_t, 32U, 0x96a5U, 16U, 11U, 19U> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::gjrand_engine<std::uint64_t, 32U, 0x96a5U, 16U, 11U, 19U> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::gjrand_engine<std::uint64_t, 32U, 0x96a5U, 16U, 11U, 19U> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::gjrand_engine<std::uint64_t, 32U, 0x96a5U, 16U, 11U, 19U> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    TEST_END;
}
