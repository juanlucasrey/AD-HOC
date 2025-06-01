#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/pcg_basic.h"

#include <PCG_engine.hpp>

auto main() -> int {

    {
        adhoc::PCG_engine rng;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = pcg32_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::PCG_engine rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::PCG_engine rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::PCG_engine rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::PCG_engine rng2;
        EXPECT_EQUAL(rng, rng2);
    }
}
