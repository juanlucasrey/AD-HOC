#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "philox_engine.hpp"

int main() {

    TEST_START;

    // std check
    {
        adhoc::philox4x32 rng1; // overload (1)
        adhoc::philox4x64 rng2; // overload (1)
        rng1.discard(10000 - 1);
        rng2.discard(10000 - 1);
        EXPECT_EQUAL(rng1(), 1955073260U);
        EXPECT_EQUAL(rng2(), 3409172418970261260U);
    }

    // check full cycle equality
    {
        adhoc::philox4x64 rng1;
        adhoc::philox4x64 rng2;
        rng2.set_counter(std::array<std::uint_fast64_t, 4UL>{
            18446744073709551615U, 18446744073709551615U, 18446744073709551615U,
            18446744073709551615U});
        for (unsigned int i = 0; i < 4; ++i) {
            rng2();
        }
        EXPECT_EQUAL(rng1, rng2);
    }

    // 32 fwd and back
    {
        adhoc::philox4x32 rng;
        adhoc::philox4x32 rng2;

        TEST_FUNC(check_fwd_and_back(rng, 1000));
        EXPECT_EQUAL(rng, rng2);

        TEST_FUNC(check_back_and_fwd(rng, 1000));
        EXPECT_EQUAL(rng, rng2);
    }

    // 64 fwd and back
    {
        adhoc::philox4x64 rng;
        adhoc::philox4x64 rng2;

        TEST_FUNC(check_fwd_and_back(rng, 1000));
        EXPECT_EQUAL(rng, rng2);

        TEST_FUNC(check_back_and_fwd(rng, 1000));
        EXPECT_EQUAL(rng, rng2);
    }

    // check 32->64
    {
        adhoc::philox4x32 rng1; // overload (1)
        adhoc::philox_engine<std::uint64_t, 32, 4, 10, 0xD2511F53, 0x9E3779B9,
                             0xCD9E8D57, 0xBB67AE85>
            rng2; // overload (1)
        TEST_FUNC(compare_rng(rng1, rng2, 100000));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    TEST_END;
}
