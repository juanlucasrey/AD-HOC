#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <lxm_engine.hpp>

#include "read_csv.hpp"

int main() {

    {
        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/lxm_vals.txt");

        auto const input =
            adhoc::readCSV<std::uint64_t>("./randomgen/lxm_state.txt");

        std::size_t j = 0;
        adhoc::lxm_engine<std::uint64_t, 64> rng(input[0], input[2], input[3],
                                                 input[4], input[5]);

        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = values_from_python[j++];
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::lxm_engine<std::uint64_t, 64> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::lxm_engine<std::uint64_t, 64> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::lxm_engine<std::uint64_t, 64> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::lxm_engine<std::uint64_t, 64> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    TEST_END;
}
