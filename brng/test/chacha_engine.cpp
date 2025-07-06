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
            adhoc::readCSV<std::uint64_t>("./randomgen/chacha_vals.txt"));

        auto const input =
            adhoc::readCSV<std::uint32_t>("./randomgen/chacha_state.txt");
        std::array<std::uint32_t, 8> key{};
        std::copy(input.begin(), input.begin() + 8, key.begin());
        adhoc::chacha20 rng(key);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::chacha rng1(1);
        adhoc::chacha20 rng2(1);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        chacha4r rng1;
        std::array<std::uint32_t, 8> key{4192101508, 3037000499, 3735928559,
                                         3735928559, 0,          0,
                                         3735928559, 3735928559};
        adhoc::chacha_engine<std::uint32_t, 4> rng2(key);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha5r rng1;
        std::array<std::uint32_t, 8> key{4192101508, 3037000499, 3735928559,
                                         3735928559, 0,          0,
                                         3735928559, 3735928559};
        adhoc::chacha_engine<std::uint32_t, 5> rng2(key);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha6r rng1;
        std::array<std::uint32_t, 8> key{4192101508, 3037000499, 3735928559,
                                         3735928559, 0,          0,
                                         3735928559, 3735928559};
        adhoc::chacha_engine<std::uint32_t, 6> rng2(key);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha8r rng1;
        std::array<std::uint32_t, 8> key{4192101508, 3037000499, 3735928559,
                                         3735928559, 0,          0,
                                         3735928559, 3735928559};
        adhoc::chacha_engine<std::uint32_t, 8> rng2(key);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        chacha20r rng1;
        std::array<std::uint32_t, 8> key{4192101508, 3037000499, 3735928559,
                                         3735928559, 0,          0,
                                         3735928559, 3735928559};
        adhoc::chacha20 rng2(key);
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
