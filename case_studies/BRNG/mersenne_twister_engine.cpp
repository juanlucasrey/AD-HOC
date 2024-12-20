#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "mersenne_twister_engine.hpp"

#include <random>
#include <vector>

int main() {
    TEST_START;

    // check against std
    {
        std::mt19937 rng1;
        adhoc::mt19937 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // check against std 64
    {
        std::mt19937_64 rng1;
        adhoc::mt19937_64 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // 32 fwd and back
    {
        adhoc::mt19937 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // 32 back and fwd
    {
        adhoc::mt19937 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // 64 fwd and back
    {
        adhoc::mt19937_64 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // 64 back and fwd
    {
        adhoc::mt19937_64 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    TEST_END;
}
