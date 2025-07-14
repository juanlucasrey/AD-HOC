#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/unimi/MWC192.c"
#include <mwc_engine.hpp>
#include <tools/mask.hpp>
#include <tools/modular_arithmetic.hpp>

auto main() -> int {

#ifndef _MSC_VER
    {
        adhoc::mcw192 rng(5678, 1234, 1);

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
        adhoc::mcw192 rng(1234, 5678, 1);
        check_fwd_and_back(rng, 1000000);
        adhoc::mcw192 rng2(1234, 5678, 1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcw192 rng(1234, 5678, 1);
        check_back_and_fwd(rng, 1000000);
        adhoc::mcw192 rng2(1234, 5678, 1);
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::mcw192>);

    TEST_END;
}
