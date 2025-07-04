#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <mersenne_twister_engine.hpp>
#include <xoshiro_engine.hpp>

#include "seed_seq_inserter.hpp"

#include "third_party/unimi/xoshiro512plus.c"

auto main() -> int {
    adhoc::mt19937_64 rng_seed;
    std::vector<std::uint64_t> inputs(8);
    std::ranges::generate(inputs, rng_seed);

    adhoc::seed_seq_inserter seq(inputs);
    adhoc::xoshiro512plus rng(seq);

    for (std::size_t i = 0; i < inputs.size(); ++i) {
        s[i] = inputs[i];
    }

    compare_rng(rng, next, 1000000);

    TEST_END;
}
