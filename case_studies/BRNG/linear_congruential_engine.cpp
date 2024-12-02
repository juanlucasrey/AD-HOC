#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "linear_congruential_engine.hpp"

#include <cstdint>
#include <random>

int main() {
    TEST_START;

    // check against std minstd_rand0
    {
        std::minstd_rand0 rng1;
        adhoc::minstd_rand0 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
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

    // check against std minstd_rand
    {
        std::minstd_rand rng1;
        adhoc::minstd_rand rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
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

    // check against std ZX81
    {
        std::linear_congruential_engine<std::uint_fast16_t, 75, 74, 0> rng1;
        adhoc::ZX81 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // ZX81 fwd and back
    {
        adhoc::ZX81 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // ZX81 back and fwd
    {
        adhoc::ZX81 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std ranqd1
    {
        std::linear_congruential_engine<std::uint_fast32_t, 1664525, 1013904223,
                                        0>
            rng1;
        adhoc::ranqd1 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // ranqd1 fwd and back
    {
        adhoc::ranqd1 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // ranqd1 back and fwd
    {
        adhoc::ranqd1 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std RANDU
    {
        std::linear_congruential_engine<std::uint_fast32_t, 65539, 0,
                                        2147483648>
            rng1;
        adhoc::RANDU rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // RANDU fwd and back
    {
        adhoc::RANDU rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // RANDU back and fwd
    {
        adhoc::RANDU rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std Borland
    {
        std::linear_congruential_engine<std::uint_fast32_t, 22695477, 1,
                                        2147483648>
            rng1;
        adhoc::Borland rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // Borland fwd and back
    {
        adhoc::Borland rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // Borland back and fwd
    {
        adhoc::Borland rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std Newlib
    {
        std::linear_congruential_engine<std::uint_fast64_t, 6364136223846793005,
                                        1, 9223372036854775808U>
            rng1;
        adhoc::Newlib rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // Newlib fwd and back
    {
        adhoc::Newlib rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // Newlib back and fwd
    {
        adhoc::Newlib rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    TEST_END;
}
