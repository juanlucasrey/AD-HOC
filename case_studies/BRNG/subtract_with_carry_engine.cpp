#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "subtract_with_carry_engine.hpp"

#include <cassert>
#include <cstdint>
#include <random>

int main() {
    TEST_START;

    // std check
    {
        adhoc::ranlux24_base gen24; // overload (1)
        adhoc::ranlux48_base gen48; // overload (1)
        gen24.discard(10000 - 1);
        gen48.discard(10000 - 1);
        assert(gen24() == 7'937'952);
        assert(gen48() == 61'839'128'582'725);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4, false>
            gen_check;
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4, false>
            gen_short;

        constexpr unsigned long long full_cycle =
            adhoc::subtract_with_carry_engine<unsigned int, 16, 2, 4,
                                              false>::period(); // 65280

        gen_short.discard(full_cycle);

        // different because of carry being initialised incorrectly!
        assert(gen_check != gen_short);

        gen_check.discard(full_cycle);
        gen_short.discard(full_cycle);

        // after a second full cycle they are in sync
        assert(gen_check == gen_short);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4, true>
            gen_check;
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4, true>
            gen_short;

        constexpr unsigned long long full_cycle =
            adhoc::subtract_with_carry_engine<unsigned int, 16, 2, 4,
                                              false>::period(); // 65280

        gen_short.discard(full_cycle);

        // carry is the same because initialisation is corrected
        assert(gen_check == gen_short);
    }

    // check against std ranlux24_base
    {
        std::ranlux24_base rng1;
        adhoc::ranlux24_base rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // check against std ranlux48_base
    {
        std::ranlux48_base rng1;
        adhoc::ranlux48_base rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // ranlux24_base fwd and back
    {
        adhoc::ranlux24_base rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // ranlux24_base back and fwd
    {
        adhoc::ranlux24_base rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // ranlux48_base fwd and back
    {
        adhoc::ranlux48_base rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // ranlux48_base back and fwd needs the correction!
    {
        adhoc::ranlux48_base rng;
        adhoc::ranlux48_base rng_check;
        rng.template operator()<false>();
        rng();
        EXPECT_NOT_EQUAL(rng, rng_check);
    }

    // ranlux48_base corrected back and fwd
    {
        adhoc::ranlux48_base_correct rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    TEST_END;
}
