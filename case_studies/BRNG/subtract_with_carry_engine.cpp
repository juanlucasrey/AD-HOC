#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "subtract_with_carry_engine.hpp"

#include <boost/random/ranlux.hpp>
#include <boost/random/subtract_with_carry.hpp>

#include <cstdint>
#include <random>

int main() {
    TEST_START;

    {
        adhoc::ranlux24_base gen24; // overload (1)
        adhoc::ranlux48_base gen48; // overload (1)
        gen24.discard(10000 - 1);
        gen48.discard(10000 - 1);
        assert(gen24() == 7'937'952);
        assert(gen48() == 61'839'128'582'725);
    }

    // check against std ranlux24_base
    {
        std::ranlux24_base rng1;
        adhoc::ranlux24_base rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // check against std ranlux24_base
    {
        std::ranlux48_base rng1;
        adhoc::ranlux48_base rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    TEST_END;
}
