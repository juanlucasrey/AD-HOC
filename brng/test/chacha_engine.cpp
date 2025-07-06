#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "third_party/imneme/chacha.hpp"
#include <chacha_engine.hpp>

#include "read_csv.hpp"
#include "seed_seq_inserter.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

int main() {

    {
        auto const values_from_python = adhoc::split_uint64_to_uint32(
            adhoc::readCSV<std::uint64_t>("./randomgen/ChaCha_vals.txt"));

        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint32_t>("./randomgen/ChaCha_state.txt"));

        adhoc::chacha20 rng(seq);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::chacha rng1(1);
        adhoc::chacha20 rng2(1, 0);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        chacha4r rng1;
        adhoc::chacha_engine<std::uint32_t, 4> rng2(0xb504f333f9de6484UL,
                                                    16045690984833335023ULL, 0,
                                                    16045690984833335023ULL);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha5r rng1;
        adhoc::chacha_engine<std::uint32_t, 5> rng2(0xb504f333f9de6484UL,
                                                    16045690984833335023ULL, 0,
                                                    16045690984833335023ULL);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha6r rng1;
        adhoc::chacha_engine<std::uint32_t, 6> rng2(0xb504f333f9de6484UL,
                                                    16045690984833335023ULL, 0,
                                                    16045690984833335023ULL);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha8r rng1;
        adhoc::chacha_engine<std::uint32_t, 8> rng2(0xb504f333f9de6484UL,
                                                    16045690984833335023ULL, 0,
                                                    16045690984833335023ULL);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha20r rng1;
        adhoc::chacha20 rng2(0xb504f333f9de6484UL, 16045690984833335023ULL, 0,
                             16045690984833335023ULL);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::chacha20 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::chacha20 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::chacha20 rng;
        check_back_and_fwd(rng, 16);
        adhoc::chacha20 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::chacha_engine<std::uint32_t, 20> rng1;
        adhoc::chacha_engine<std::uint64_t, 20> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::chacha20>);

    TEST_END;
}
