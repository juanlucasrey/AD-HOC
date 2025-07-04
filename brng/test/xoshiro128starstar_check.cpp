#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <mersenne_twister_engine.hpp>
#include <xoshiro_engine.hpp>

#include "seed_seq_inserter.hpp"

#include "third_party/unimi/xoshiro128starstar.c"

auto main() -> int {
    adhoc::mt19937 rng_seed;
    std::vector<std::uint32_t> inputs(4);
    std::ranges::generate(inputs, rng_seed);

    adhoc::seed_seq_inserter seq(inputs);
    adhoc::xoshiro128starstar rng(seq);

    for (std::size_t i = 0; i < inputs.size(); ++i) {
        s[i] = inputs[i];
    }

    compare_rng(rng, next, 1000000);

    TEST_END;
}
