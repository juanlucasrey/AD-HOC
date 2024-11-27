#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "linear_congruential_engine.hpp"

#include <random>
#include <vector>

int main() {
    TEST_START;

    // check against std minstd_rand0
    {
        std::minstd_rand0 rng1;
        adhoc::minstd_rand0 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // check against std minstd_rand
    {
        std::minstd_rand rng1;
        adhoc::minstd_rand rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // minstd_rand0 fwd and back
    {
        adhoc::minstd_rand0 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // minstd_rand0 back and fwd
    {
        adhoc::minstd_rand0 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // minstd_rand fwd and back
    {
        adhoc::minstd_rand rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // minstd_rand back and fwd
    {
        adhoc::minstd_rand rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    TEST_END;
}
