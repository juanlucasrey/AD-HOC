#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

extern "C" {
#include "third_party/umontreal/MRG32k3a.h"
}

#include "seed_seq_inserter.hpp"
#include <MRG32k3a.hpp>
#include <distribution/uniform_distribution.hpp>

#include <cstdint>

auto main() -> int {
    {
        // x, x, 0, x, x, 0 precedent
        std::vector<std::uint32_t> seqval{4243840401, 3757362166, 89818307,
                                          3482899528, 2594243115, 3560662590};
        adhoc::seed_seq_inserter seq(seqval);
        adhoc::open<adhoc::MRG32k3a<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {4243840401, 3757362166, 89818307,
                                           3482899528, 2594243115, 3560662590};
        InitMRG32k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // x, x, x, x, x, 0 precedent
        std::vector<std::uint32_t> seqval{1806066723, 2321277356, 1826996185,
                                          2720448313, 1027871599, 1848796065};
        adhoc::seed_seq_inserter seq(seqval);
        adhoc::open<adhoc::MRG32k3a<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {1806066723, 2321277356, 1826996185,
                                           2720448313, 1027871599, 1848796065};
        InitMRG32k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // x, x, 0, x, x, x precedent
        std::vector<std::uint32_t> seqval{2612542496, 3458909419, 1729080824,
                                          1842728981, 3072837856, 2942928930};
        adhoc::seed_seq_inserter seq(seqval);
        adhoc::open<adhoc::MRG32k3a<std::uint32_t, true>> rng(seq);

        std::vector<unsigned long> init = {2612542496, 3458909419, 1729080824,
                                           1842728981, 3072837856, 2942928930};
        InitMRG32k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::open<adhoc::MRG32k3a<std::uint32_t, true>> rng;
        std::vector<unsigned long> init = {12345, 12345, 12345,
                                           12345, 12345, 12345};
        InitMRG32k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG32k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::MRG32k3a<std::uint32_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG32k3a<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k3a<std::uint32_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG32k3a<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k3a<std::uint64_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG32k3a<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k3a<std::uint64_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG32k3a<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG32k3a<std::uint32_t> rng1;
        adhoc::MRG32k3a<std::uint64_t> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::MRG32k3a<std::uint32_t>>);

    TEST_END;
}
