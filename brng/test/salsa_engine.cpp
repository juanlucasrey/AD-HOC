#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/imneme/chacha.hpp"
#include <salsa_engine.hpp>

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

int main() {

    {
        PractRand::RNGs::LightWeight::salsa rng1(0xb504f333f9de6484UL);
        adhoc::salsa20 rng2;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }
    {
        PractRand::RNGs::LightWeight::salsa rng1(1);
        adhoc::salsa20 rng2(1);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::salsa20 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::salsa20 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::salsa_engine<std::uint32_t, 20> rng1;
        adhoc::salsa_engine<std::uint64_t, 20> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    TEST_END;
}
