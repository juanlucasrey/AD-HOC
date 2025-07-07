#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <squares_engine.hpp>

#include "read_csv.hpp"

int main() {

    {
        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/squares_vals.txt");

        auto const input =
            adhoc::readCSV<std::uint64_t>("./randomgen/squares_state.txt");
        std::uint64_t key = input[0];
        std::size_t j = 0;
        adhoc::squares64 rng(key);
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        auto const values_from_python = adhoc::split_uint64_to_uint32(
            adhoc::readCSV<std::uint64_t>("./randomgen/squares_32_vals.txt"));

        auto const input =
            adhoc::readCSV<std::uint64_t>("./randomgen/squares_32_state.txt");
        std::uint64_t key = input[0];
        adhoc::squares32 rng(key);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::squares64 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::squares64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::squares64 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::squares64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::squares32 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::squares32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::squares32 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::squares32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::squares64>);
    static_assert(std::bidirectional_iterator<adhoc::squares32>);
    TEST_END;
}
