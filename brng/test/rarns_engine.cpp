#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <rarns_engine.hpp>

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

auto main() -> int {

    {
        PractRand::RNGs::LightWeight::rarns64 rng1(1);
        adhoc::rarns64<std::uint64_t> rng2(1U, 1U, 1U, 8);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw64();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::rarns32 rng1(1);
        adhoc::rarns32<std::uint32_t> rng2(1U, 0U);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::rarns16 rng1(1);
        adhoc::rarns16<std::uint16_t> rng2(1U, 0U, 0U, 0);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw16();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::rarns64<std::uint64_t> rng(1U, 1U, 1U, 8);
        check_fwd_and_back(rng, 1000000);
        adhoc::rarns64<std::uint64_t> rng2(1U, 1U, 1U, 8);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::rarns64<std::uint64_t> rng(1U, 1U, 1U, 8);
        check_back_and_fwd(rng, 1000000);
        adhoc::rarns64<std::uint64_t> rng2(1U, 1U, 1U, 8);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::rarns32<std::uint32_t> rng(1U, 0U);
        check_fwd_and_back(rng, 1000000);
        adhoc::rarns32<std::uint32_t> rng2(1U, 0U);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::rarns32<std::uint32_t> rng(1U, 0U);
        check_back_and_fwd(rng, 1000000);
        adhoc::rarns32<std::uint32_t> rng2(1U, 0U);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::rarns16<std::uint16_t> rng(1U, 0U, 0U, 0);
        check_fwd_and_back(rng, 1000000);
        adhoc::rarns16<std::uint16_t> rng2(1U, 0U, 0U, 0);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::rarns_engine<std::uint32_t, 32, 10, 3, 5, 13> rng1;
        adhoc::rarns_engine<std::uint64_t, 32, 10, 3, 5, 13> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    TEST_END;
}
