#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "../../include/combinatorics/pow.hpp"
#include "subtract_with_carry_engine.hpp"

#include <cassert>
#include <cstdint>
#include <random>

template <class RNG>
constexpr auto period(RNG const & /* rng */) -> unsigned long long {
    constexpr auto w = RNG::word_size;
    constexpr auto r = RNG::long_lag;
    constexpr auto s = RNG::short_lag;
    return adhoc::detail::pow<r>(w) - adhoc::detail::pow<s>(w);
}

class initseq : public std::seed_seq {
  public:
    template <class RandomIt> void generate(RandomIt first, RandomIt last) {
        if (first == last)
            return;

        std::size_t i = 0;
        for (; first != last; ++first) {
            *first = vals[i];
            i++;
        }
    }
    std::vector<std::uint_least32_t> vals;
};

auto init_test(std::vector<std::uint_least32_t> &&init) -> int {
    TEST_START;

    initseq seed;
    seed.vals = init;

    adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng1(seed);
    adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng2(seed);
    constexpr unsigned long long full_cycle = period(rng1); // 65280
    rng2.discard(full_cycle);
    EXPECT_EQUAL(rng1, rng2);

    TEST_FUNC(check_fwd_and_back(rng2, full_cycle));
    TEST_FUNC(check_back_and_fwd(rng2, full_cycle));

    TEST_END;
}

int main() {
    TEST_START;

    // std check
    {
        adhoc::ranlux24_base gen24; // overload (1)
        adhoc::ranlux48_base gen48; // overload (1)
        gen24.discard(10000 - 1);
        gen48.discard(10000 - 1);
        EXPECT_EQUAL(gen24(), 7'937'952);
        EXPECT_EQUAL(gen48(), 61'839'128'582'725);
    }

    // to be corrected on std!!
    {
        std::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng1;
        std::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng2;

        constexpr unsigned long long full_cycle = period(rng1); // 65280

        rng2.discard(full_cycle);

        // std flaw!! this should be equal
        EXPECT_NOT_EQUAL(rng1, rng2);

        // but they give the same values!!!
        TEST_FUNC(compare_rng(rng1, rng2, full_cycle));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // default
    {
        initseq seed;
        seed.vals = std::vector<std::uint_least32_t>{63026, 2533, 52484, 29490};
        std::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng1(
            seed);
        std::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng2;
        EXPECT_EQUAL(rng1, rng2);
    }

    // underflow and overflow inits
    {
        TEST_FUNC(init_test(
            std::vector<std::uint_least32_t>{63026, 2533, 52484, 29490}));

        // carry > carry_std with underflow
        TEST_FUNC(
            init_test(std::vector<std::uint_least32_t>{0, 2533, 52484, 29490}));

        // carry > carry_std with 2 underflows
        TEST_FUNC(
            init_test(std::vector<std::uint_least32_t>{0, 0, 52484, 29490}));

        // carry > carry_std with 3 underflows
        TEST_FUNC(init_test(std::vector<std::uint_least32_t>{0, 0, 0, 29490}));

        // carry > carry_std with 4 underflows? all 0 makes carry = 1
        TEST_FUNC(init_test(std::vector<std::uint_least32_t>{0, 0, 0, 0}));

        // carry < carry_std
        TEST_FUNC(
            init_test(std::vector<std::uint_least32_t>{63026, 2533, 52484, 0}));

        // carry < carry_std one overflow
        TEST_FUNC(
            init_test(std::vector<std::uint_least32_t>{65535, 2533, 52484, 0}));

        // carry < carry_std two overflow
        TEST_FUNC(init_test(
            std::vector<std::uint_least32_t>{65535, 65535, 52484, 0}));

        // carry < carry_std three overflow
        TEST_FUNC(init_test(
            std::vector<std::uint_least32_t>{65535, 65535, 65535, 0}));

        // carry < carry_std one overflow
        TEST_FUNC(
            init_test(std::vector<std::uint_least32_t>{65535, 2533, 52484, 1}));

        // carry < carry_std two overflow
        TEST_FUNC(init_test(
            std::vector<std::uint_least32_t>{65535, 65535, 52484, 1}));

        // carry < carry_std three overflow
        TEST_FUNC(init_test(
            std::vector<std::uint_least32_t>{65535, 65535, 65535, 1}));
    }

    // subtract_with_carry_engine with small period
    {
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4>
            gen_check;
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4>
            gen_short;

        constexpr unsigned long long full_cycle = period(gen_check); // 65280

        gen_short.discard(full_cycle);
        assert(gen_check == gen_short);
    }

    // check against std ranlux24_base
    {
        std::ranlux24_base rng1;
        adhoc::ranlux24_base rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
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

    // check against std ranlux48_base
    {
        std::ranlux48_base rng1;
        adhoc::ranlux48_base rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
        TEST_FUNC(compare_rng_limits(rng1, rng2));
    }

    // ranlux48_base fwd and back
    {
        adhoc::ranlux48_base rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // ranlux48_base corrected back and fwd
    {
        adhoc::ranlux48_base rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    TEST_END;
}
