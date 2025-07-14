#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/unimi/MWC128.c"
#include <mwc_engine.hpp>
#include <tools/mask.hpp>
#include <tools/modular_arithmetic.hpp>

auto main() -> int {

    {
        adhoc::mcw128 rng(1234, 1);

        x = 1234;
        c = 1;
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = next();
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::mcw128 rng(1234, 1);
        check_fwd_and_back(rng, 1000000);
        adhoc::mcw128 rng2(1234, 1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcw128 rng(1234, 1);
        check_back_and_fwd(rng, 1000000);
        adhoc::mcw128 rng2(1234, 1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcw128 rng(1234, 7045691847467367358);
        ++rng;
        ++rng;
        // the rng now holds a state of c=0, so I don't know why there's such a
        // fuss about exluding the case c=0 from the initial state..

        // to prove that c=0 we show that we have x_{i+1} = x_i * A (implying c
        // = 0)
        auto val1 = *rng;
        ++rng;
        auto val2 = *rng;
        constexpr uint64_t MWC_A = 0xffebb71d94fcdaf9;
        auto val2b = val1 * MWC_A;
        EXPECT_EQUAL(val2, val2b);
    }

    static_assert(std::bidirectional_iterator<adhoc::mcw128>);

    TEST_END;
}
