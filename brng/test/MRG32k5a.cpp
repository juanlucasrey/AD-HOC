#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

extern "C" {
#include "third_party/umontreal/MRG32k5a.h"
}

#include "seed_seq_inserter.hpp"
#include <MRG32k5a.hpp>
#include <distribution/uniform_distribution.hpp>

#include <cstdint>

auto main() -> int {

    {
        // p2 == 0 case
        std::vector<std::uint32_t> seqval{
            1958623832, 2866876007, 150144429,  478980109, 740639212,
            3955543965, 2550687073, 3509537915, 119502490, 4039470558};
        adhoc::seed_seq_inserter seq(seqval);
        adhoc::open<adhoc::MRG32k5a<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {
            1958623832, 2866876007, 150144429,  478980109, 740639212,
            3955543965, 2550687073, 3509537915, 119502490, 4039470558};
        InitMRG32k5a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k5a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // p1 == 0 case
        std::vector<std::uint32_t> seqval{
            1958623832, 2866876007, 150144429, 478980109,  740639212,
            2313483013, 1786764733, 713760824, 4050934279, 1671982563};
        adhoc::seed_seq_inserter seq(seqval);
        adhoc::open<adhoc::MRG32k5a<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {
            1958623832, 2866876007, 150144429, 478980109,  740639212,
            2313483013, 1786764733, 713760824, 4050934279, 1671982563};
        InitMRG32k5a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k5a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // p1 == p2 case
        std::vector<std::uint32_t> seqval{
            2027557476, 2763277694, 1990012522, 2412316413, 2112525205,
            2313483013, 1786764733, 713760824,  4050934279, 1671982563};
        adhoc::seed_seq_inserter seq(seqval);
        adhoc::open<adhoc::MRG32k5a<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {
            2027557476, 2763277694, 1990012522, 2412316413, 2112525205,
            2313483013, 1786764733, 713760824,  4050934279, 1671982563};
        InitMRG32k5a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k5a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::open<adhoc::MRG32k5a<std::uint32_t, true>> rng;
        std::vector<unsigned long> init = {12345, 12345, 12345, 12345, 12345,
                                           12345, 12345, 12345, 12345, 12345};
        InitMRG32k5a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k5a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::MRG32k5a<std::uint32_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG32k5a<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k5a<std::uint32_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG32k5a<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k5a<std::uint64_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG32k5a<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k5a<std::uint64_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG32k5a<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k5a<std::uint32_t> rng1;
        adhoc::MRG32k5a<std::uint64_t> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::MRG32k5a<std::uint32_t>>);

    TEST_END;
}
