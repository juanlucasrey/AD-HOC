#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "philox_engine.hpp"

// #include <counter_based_engine.hpp>
#include <alternative_version/philox_engine.hpp>

#include <random>

int main() {
    TEST_START;

    // check against std
    {
        adhoc::philox4x32 rng1; // overload (1)
        adhoc::philox4x64 rng2; // overload (1)
        rng1.discard(10000 - 1);
        rng2.discard(10000 - 1);
        EXPECT_EQUAL(rng1(), 1955073260U);
        EXPECT_EQUAL(rng2(), 3409172418970261260U);
    }

    // literal std description does not work!
    {
        adhoc::philox4x32 rng1; // overload (1)
        adhoc::philox4x64 rng2; // overload (1)
        rng1.discard<0>(10000 - 1);
        rng2.discard<0>(10000 - 1);
        EXPECT_NOT_EQUAL(rng1(), 1955073260U);
        EXPECT_NOT_EQUAL(rng2(), 3409172418970261260U);
    }

    // check against ref impl
    {
        philox4x32 rng1;        // overload (1)
        adhoc::philox4x32 rng2; // overload (1)
        TEST_FUNC(compare_rng(rng1, rng2, 100000));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
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

    // check against ref impl 32->64
    {
        philox_engine<std::uint64_t, 32, 4, 10, 1, 0xD2511F53, 0x9E3779B9,
                      0xCD9E8D57, 0xBB67AE85>
            rng1; // overload (1)
        adhoc::philox_engine<std::uint64_t, 32, 4, 10, 0xD2511F53, 0x9E3779B9,
                             0xCD9E8D57, 0xBB67AE85>
            rng2; // overload (1)
        TEST_FUNC(compare_rng(rng1, rng2, 100000));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // check against ref impl
    {
        philox4x64 rng1;        // overload (1)
        adhoc::philox4x64 rng2; // overload (1)
        TEST_FUNC(compare_rng(rng1, rng2, 100000));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // check against ref impl
    {
        philox_engine<std::uint_fast32_t, 32, 2, 10, 1, 0xD2511F53,
                      0x9E3779B9>
            rng1; // overload (1)
        adhoc::philox_engine<std::uint_fast32_t, 32, 2, 10, 0xD2511F53,
                             0x9E3779B9>
            rng2; // overload (1)
        TEST_FUNC(compare_rng(rng1, rng2, 100000));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // check against ref impl
    {
        philox_engine<std::uint_fast64_t, 64, 2, 10, 1, 0xD2E7470EE14C6C93,
                      0x9E3779B97F4A7C15>
            rng1; // overload (1)
        adhoc::philox_engine<std::uint_fast64_t, 64, 2, 10, 0xD2E7470EE14C6C93,
                             0x9E3779B97F4A7C15>
            rng2; // overload (1)
        TEST_FUNC(compare_rng(rng1, rng2, 100000));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    TEST_END;
}
