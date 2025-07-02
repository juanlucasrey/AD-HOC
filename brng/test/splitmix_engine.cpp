#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/imneme/splitmix.hpp"
#include <splitmix_engine.hpp>

int main() {

    {
        splitmix64 rng1;
        adhoc::splitmix64 rng2;

        rng1.advance(10004);
        rng2.discard(10004);

        auto rng1b = rng1.split();
        auto rng2b = rng2.split();
        compare_rng(rng1b, rng2b, 10000);
        compare_rng_limits(rng1b, rng2b);
    }

    {
        splitmix32 rng1;
        adhoc::splitmix32 rng2;
        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::splitmix32 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::splitmix32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::splitmix32 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::splitmix32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::splitmix32 rng1;
        adhoc::splitmix_engine<
            std::uint64_t, 32, 0xff51afd7ed558ccdUL, 0xc4ceb9fe1a85ec53UL, 33,
            33, 33, 0x62a9d9ed799705f5UL, 0xcb24d0a5c88c35b3UL, 33, 28, 0>
            rng2;
        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        splitmix64 rng1;
        adhoc::splitmix64 rng2;
        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::splitmix64 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::splitmix64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::splitmix64 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::splitmix64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::splitmix64 rng1;
        adhoc::splitmix64 rng2;
        for (std::size_t i = 0; i < 1000000; ++i) {
            rng2();
        }
        auto diff = rng2 - rng1;
        EXPECT_EQUAL(diff, 1000000);
    }

    static_assert(std::bidirectional_iterator<adhoc::splitmix64>);
    static_assert(std::bidirectional_iterator<adhoc::splitmix32>);

    TEST_END;
}
