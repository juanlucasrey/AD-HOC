#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/imneme/chacha.hpp"
#include <chacha_engine.hpp>

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

int main() {

    {
        PractRand::RNGs::LightWeight::chacha rng1(1);
        adhoc::chacha_engine<std::uint32_t, 32, 20, true> rng2(1, 0);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        chacha4r rng1;
        adhoc::chacha4r rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha5r rng1;
        adhoc::chacha5r rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha6r rng1;
        adhoc::chacha6r rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha8r rng1;
        adhoc::chacha8r rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha20r rng1;
        adhoc::chacha20r rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::chacha4r rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::chacha4r rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::chacha4r rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::chacha4r rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    TEST_END;
}
