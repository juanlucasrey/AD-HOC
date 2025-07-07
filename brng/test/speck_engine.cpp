#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "seed_seq_inserter.hpp"
#include <speck_engine.hpp>

#include "read_csv.hpp"

int main() {

    {
        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/speck_vals.txt");

        auto const input =
            adhoc::readCSV<std::uint64_t>("./randomgen/speck_state.txt");

        std::vector<std::uint64_t> input_2(34);
        for (std::size_t i = 0; i < 34; i++) {
            input_2[i] = input[i * 2];
        }

        adhoc::seed_seq_inserter seq(input_2);

        adhoc::speck34 rng(seq);

        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = values_from_python[j++];
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    TEST_END;
}
