#include "../../test_simple/test_simple_include.hpp"

#include "RNGtest.hpp"

#include <distribution/normal_distribution.hpp>
#include <mersenne_twister_engine.hpp>

#include <cstdint>
#include <map>
#include <vector>

auto main() -> int {
    {
        constexpr std::uint32_t min = 0;
        constexpr std::uint32_t max = std::numeric_limits<std::uint32_t>::max();
        constexpr double largest_normal = 6.2302601380255185;

        {
            adhoc::standard_normal<RNGtest<std::uint32_t, min, max>> gen(max);
            auto val = gen();
            EXPECT_EQUAL(val, largest_normal);
        }

        {
            adhoc::standard_normal<RNGtest<std::uint32_t, min, max>> gen(min);
            auto val = gen();
            EXPECT_EQUAL(val, -largest_normal);
        }
    }

    {
        constexpr std::uint64_t min = 0;
        constexpr std::uint64_t max = std::numeric_limits<std::uint64_t>::max();

        constexpr double largest_normal = 9.0801551248736114;

        {
            adhoc::standard_normal<RNGtest<std::uint64_t, min, max>> gen(max);
            auto val = gen();
            EXPECT_EQUAL(val, largest_normal);
        }

        {
            adhoc::standard_normal<RNGtest<std::uint64_t, min, max>> gen(min);
            auto val = gen();
            EXPECT_EQUAL(val, -largest_normal);
        }
    }

#ifndef _MSC_VER
    {
        constexpr __uint128_t min = 0;
        constexpr __uint128_t max = std::numeric_limits<__uint128_t>::max();

        constexpr double largest_normal = 13.055946840098045;

        {
            adhoc::standard_normal<RNGtest<__uint128_t, min, max>> gen(max);
            auto val = gen();
            EXPECT_EQUAL(val, largest_normal);
        }

        {
            adhoc::standard_normal<RNGtest<__uint128_t, min, max>> gen(min);
            auto val = gen();
            EXPECT_EQUAL(val, -largest_normal);
        }
    }
#endif

    {
        adhoc::standard_normal<adhoc::mt19937_64> gen;

        // Draw a sample from the normal distribution and round it to an
        // integer.
        auto random_int = [&gen] { return std::lround(gen()); };

        std::map<long, unsigned> histogram{};
        for (auto n{10000}; n; --n) {
            ++histogram[random_int()];
        }

        std::vector<unsigned> result{3, 63, 598, 2448, 3780, 2441, 593, 70, 4};
        std::size_t i = 0;
        for (const auto [k, v] : histogram) {
            auto res = result[i++];
            EXPECT_EQUAL(v, res);
        }
    }

    TEST_END;
}
