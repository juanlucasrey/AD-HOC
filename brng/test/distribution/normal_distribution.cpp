#include "../../test_simple/test_simple_include.hpp"

#include "RNGtest.hpp"

#include <distribution/normal_distribution.hpp>
#include <mersenne_twister_engine.hpp>

#include <cstdint>
#include <map>
#include <vector>

auto main() -> int {
    {
        constexpr std::uint64_t min64 = 0;
        constexpr std::uint64_t max64 =
            std::numeric_limits<std::uint64_t>::max();

        constexpr double largest_normal = 9.0801551248736114;

        {
            adhoc::standard_normal<RNGtest<std::uint64_t, min64, max64>> gen(
                max64);
            auto val = gen();
            EXPECT_EQUAL(val, largest_normal);
        }

        {
            adhoc::standard_normal<RNGtest<std::uint64_t, min64, max64>> gen(
                min64);
            auto val = gen();
            EXPECT_EQUAL(val, -largest_normal);
        }
    }

    {
        constexpr std::uint32_t min32 = 0;
        constexpr std::uint32_t max32 =
            std::numeric_limits<std::uint32_t>::max();
        constexpr double largest_normal = 6.2302601380255185;

        {
            adhoc::standard_normal<RNGtest<std::uint32_t, min32, max32>> gen(
                max32);
            auto val = gen();
            EXPECT_EQUAL(val, largest_normal);
        }

        {
            adhoc::standard_normal<RNGtest<std::uint32_t, min32, max32>> gen(
                min32);
            auto val = gen();
            EXPECT_EQUAL(val, -largest_normal);
        }
    }

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
