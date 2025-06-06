#include "../../../test_simple/test_simple_include.hpp"

#include "../include/tools/bit_find.hpp"

#include <cstdint>
#include <random>

auto max64(int w) -> std::uint64_t {
    uint64_t ach = static_cast<uint64_t>(1U)
                   << static_cast<std::size_t>(w - 1U);

    return ach | (ach - static_cast<uint64_t>(1U));

    // equivalent to the following but shift does not overflow
    // return (static_cast<result_type>(1U) << w) - 1U;
}

auto max32(int w) -> std::uint32_t {
    uint32_t ach = std::uint32_t{1U} << static_cast<std::size_t>(w - 1U);

    return ach | (ach - std::uint32_t{1U});

    // equivalent to the following but shift does not overflow
    // return (static_cast<result_type>(1U) << w) - 1U;
}

int main() {
    {
        for (std::uint32_t i = 0; i < 1000; i++) {
            auto res1 = adhoc::findLeastSignificantZero<false>(i);
            auto res2 = adhoc::findLeastSignificantZero<true>(i);
            EXPECT_EQUAL(res1, res2);
        }
    }

    {
        for (std::size_t i = 0; i < 32; i++) {
            std::uint32_t val = max32(static_cast<int>(i));
            auto res1 = adhoc::findLeastSignificantZero<false>(val);
            auto res2 = adhoc::findLeastSignificantZero<true>(val);
            EXPECT_EQUAL(res1, res2);
        }
    }

    {
        std::mt19937_64 gen(123);
        std::uniform_int_distribution<uint64_t> dis;

        for (std::size_t i = 0; i < 1000; i++) {
            std::uint64_t val = dis(gen);
            auto res1 = adhoc::findLeastSignificantZero<false>(val);
            auto res2 = adhoc::findLeastSignificantZero<true>(val);
            EXPECT_EQUAL(res1, res2);
        }
    }

    {
        for (std::size_t i = 0; i < 64; i++) {
            std::uint64_t val = max64(static_cast<int>(i));
            auto res1 = adhoc::findLeastSignificantZero<false>(val);
            auto res2 = adhoc::findLeastSignificantZero<true>(val);
            EXPECT_EQUAL(res1, res2);
        }
    }

    TEST_END;
}
