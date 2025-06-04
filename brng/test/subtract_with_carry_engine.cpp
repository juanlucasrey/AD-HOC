#include "../../test_simple/test_simple_include.hpp"
#include "seed_seq.hpp"
#include "test_tools_rng.hpp"

#include "../include/subtract_with_carry_engine.hpp"

#include <cassert>
#include <chrono>
#include <cstdint>
#include <random>

namespace {

template <std::size_t N, class T> inline constexpr auto pow(T base) -> T {
    if constexpr (N == 1) {
        return base;
    } else if constexpr (N % 2 == 0) {
        T power = pow<N / 2>(base);
        return power * power;

    } else {
        T power = pow<N / 2>(base);
        return base * power * power;
    }
};

template <class RNG>
constexpr auto period(RNG const & /* rng */) -> unsigned long long {
    constexpr auto w = RNG::word_size;
    constexpr auto r = RNG::long_lag;
    constexpr auto s = RNG::short_lag;
    return pow<r>(w) - pow<s>(w);
}

auto init_test(std::vector<std::uint_least32_t> &&init) -> int {
    adhoc::seed_seq<std::uint_least32_t> seq;
    seq.vals = init;

    adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng1(seq);
    adhoc::subtract_with_carry_engine<std::uint_fast32_t, 16, 2, 4> rng2(seq);
    constexpr unsigned long long full_cycle = period(rng1); // 65280
    rng2.discard(full_cycle);
    EXPECT_EQUAL(rng1, rng2);

    check_fwd_and_back(rng2, full_cycle);
    check_back_and_fwd(rng2, full_cycle);

    TEST_END;
}

} // namespace

int main() {
    // std check
    {
        adhoc::ranlux24_base gen24; // overload (1)
        adhoc::ranlux48_base gen48; // overload (1)
        gen24.discard(10000 - 1);
        gen48.discard(10000 - 1);
        EXPECT_EQUAL(gen24(), 7'937'952);
        EXPECT_EQUAL(gen48(), 61'839'128'582'725);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 32, 10, 24> rng1;
        std::subtract_with_carry_engine<std::uint_fast32_t, 32, 10, 24> rng2;
        compare_rng(rng1, rng2, 1000000);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast64_t, 64, 5, 12> rng1;
        std::subtract_with_carry_engine<std::uint_fast64_t, 64, 5, 12> rng2;
        compare_rng(rng1, rng2, 1000000);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 32, 10, 24> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 32, 10, 24> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 32, 10, 24> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::subtract_with_carry_engine<std::uint_fast32_t, 32, 10, 24> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast64_t, 64, 5, 12> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::subtract_with_carry_engine<std::uint_fast64_t, 64, 5, 12> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::subtract_with_carry_engine<std::uint_fast64_t, 64, 5, 12> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::subtract_with_carry_engine<std::uint_fast64_t, 64, 5, 12> rng2;
        EXPECT_EQUAL(rng, rng2);
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
        compare_rng(rng1, rng2, full_cycle);
        compare_rng_limits(rng1, rng2);
    }

    // default
    {
        adhoc::seed_seq<std::uint_least32_t> seed;
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
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        std::ranlux48_base rng1;
        adhoc::ranlux48_base rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    for (std::size_t i = 0; i < 3; i++) {

        {
            adhoc::ranlux24_base rng(static_cast<std::uint_fast32_t>(i));
            check_fwd_and_back(rng, 1000000);
            adhoc::ranlux24_base rng2(static_cast<std::uint_fast32_t>(i));
            EXPECT_EQUAL(rng, rng2);
        }

        {
            adhoc::ranlux24_base rng(static_cast<std::uint_fast32_t>(i));
            check_back_and_fwd(rng, 1000000);
            adhoc::ranlux24_base rng2(static_cast<std::uint_fast32_t>(i));
            EXPECT_EQUAL(rng, rng2);
        }

        {
            adhoc::ranlux48_base rng(static_cast<std::uint_fast64_t>(i));
            check_fwd_and_back(rng, 1000000);
            adhoc::ranlux48_base rng2(static_cast<std::uint_fast64_t>(i));
            EXPECT_EQUAL(rng, rng2);
        }

        {
            adhoc::ranlux48_base rng(static_cast<std::uint_fast64_t>(i));
            check_back_and_fwd(rng, 1000000);
            adhoc::ranlux48_base rng2(static_cast<std::uint_fast64_t>(i));
            EXPECT_EQUAL(rng, rng2);
        }
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::ranlux48_base rng;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new" << std::endl;
            std::cout << time << std::endl;
        }

        std::uint64_t res2 = 0;
        std::ranlux48_base rng2;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res2 += rng2();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "old" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res2);

        std::uint64_t res1b = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1b += rng.operator()<false>();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res1b);
    }

    TEST_END;
}
