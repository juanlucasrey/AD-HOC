#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/imneme/lehmer.hpp"
#include <lehmer_engine.hpp>
#include <linear_congruential_engine.hpp>

auto main() -> int {

    {
        mcg128 rng1;
        adhoc::mcg128 rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        mcg128_fast rng1;
        adhoc::mcg128_fast rng2;

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::mcg128 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::mcg128 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcg128 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::mcg128 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcg128_fast rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::mcg128_fast rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::mcg128_fast rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::mcg128_fast rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::lehmer_engine<std::uint64_t, 32, std::uint64_t, 32, 0x42e4d8b5UL>
            rng1;
        adhoc::lehmer_engine<std::uint32_t, 32, std::uint32_t, 32, 0x42e4d8b5UL>
            rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::lehmer_engine<std::uint16_t, 16, std::uint32_t, 32, 0x42e4d8b5UL>
            rng1;
        adhoc::lehmer_engine<std::uint32_t, 16, std::uint64_t, 32, 0x42e4d8b5UL>
            rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::lehmer_engine<std::uint64_t, 64, adhoc::uint128, 128,
                             0xda942042e4dd58b5ULL>
            rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::lehmer_engine<std::uint64_t, 64, adhoc::uint128, 128,
                             0xda942042e4dd58b5ULL>
            rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::mcg128>);

    TEST_END;
}
