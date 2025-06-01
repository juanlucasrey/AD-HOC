#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

extern "C" {
#include "external/MRG31k3p.h"
}

#include "seed_seq.hpp"
#include <MRG31k3p.hpp>
#include <distribution/uniform_distribution.hpp>

#include <cstdint>

auto main() -> int {
    {
        // 0, 0
        adhoc::seed_seq<std::uint_fast32_t> seq;
        seq.vals = {1618958054, 267881243,  156255863,
                    1730679216, 1079478356, 1334495809};
        adhoc::open<adhoc::MRG31k3p<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {1618958054, 267881243,  156255863,
                                           1730679216, 1079478356, 1334495809};
        InitMRG31k3p(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG31k3p();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // 0, x
        adhoc::seed_seq<std::uint_fast32_t> seq;
        seq.vals = {856965515, 1525874275, 784371780,
                    313650421, 1918098087, 748752505};
        adhoc::open<adhoc::MRG31k3p<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {856965515, 1525874275, 784371780,
                                           313650421, 1918098087, 748752505};
        InitMRG31k3p(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG31k3p();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // x, 0
        adhoc::seed_seq<std::uint_fast32_t> seq;
        seq.vals = {122304568,  1137184201, 506722251,
                    1815112525, 291659175,  431495294};
        adhoc::open<adhoc::MRG31k3p<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {122304568,  1137184201, 506722251,
                                           1815112525, 291659175,  431495294};
        InitMRG31k3p(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG31k3p();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::open<adhoc::MRG31k3p<std::uint32_t, true>> rng;
        std::vector<unsigned long> init = {12345, 12345, 12345,
                                           12345, 12345, 12345};
        InitMRG31k3p(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG31k3p();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::MRG31k3p<std::uint32_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG31k3p<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint32_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG31k3p<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint64_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG31k3p<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint64_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG31k3p<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint32_t> rng1;
        adhoc::MRG31k3p<std::uint64_t> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    TEST_END;
}
