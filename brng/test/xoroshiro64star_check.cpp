#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <mersenne_twister_engine.hpp>
#include <xoroshiro_engine.hpp>

#include "seed_seq_inserter.hpp"

#include "third_party/unimi/xoroshiro64star.c"

#include <algorithm>

auto main() -> int {
    adhoc::mt19937 rng_seed;
    std::vector<std::uint32_t> inputs(2);
    std::ranges::generate(inputs, rng_seed);

    adhoc::seed_seq_inserter seq(inputs);
    adhoc::xoroshiro64star rng(seq);

    for (std::size_t i = 0; i < inputs.size(); ++i) {
        s[i] = inputs[i];
    }

    compare_rng(rng, next, 1000000);

    TEST_END;
}
