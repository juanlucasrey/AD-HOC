#include "../../test_simple/test_simple_include.hpp"

#include "RNGtest.hpp"

#include <distribution/centered_canonical_distribution.hpp>
#include <distribution/exponential_distribution.hpp>

#include <cstdint>
#include <iostream>

auto main() -> int {

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    constexpr double small_val = 5.4210108624275222e-20;

    // 0, max
    {
        constexpr std::uint64_t min = 0;
        constexpr std::uint64_t max = std::numeric_limits<std::uint64_t>::max();

        // open, double, max
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // open, float, max
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>,
                                           float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // open, double, min
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val);
        }

        // open, float, min
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>,
                                           float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val);
        }

        // closed, double, max
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, float, max
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>,
                                             float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, double, min
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // closed, float, min
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>,
                                             float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // default, double, max
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // default, float, max
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>,
                                              float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // default, double, min
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // default, float, min
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>,
                                              float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }
    }

    // open, 1, max
    {
        constexpr std::uint64_t min = 2;
        constexpr std::uint64_t max = std::numeric_limits<std::uint64_t>::max();

        // open, double, max
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // open, float, max
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>,
                                           float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // open, double, min
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val);
        }

        // open, float, min
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>,
                                           float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val);
        }

        // closed, double, max
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, float, max
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>,
                                             float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, double, min
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // closed, float, min
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>,
                                             float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // default, double, max
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // default, float, max
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>,
                                              float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // default, double, min
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // default, float, min
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>,
                                              float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }
    }

    // open, 2, max
    {
        constexpr std::uint64_t min = 1;
        constexpr std::uint64_t max = std::numeric_limits<std::uint64_t>::max();

        // open, double, max
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // open, float, max
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>,
                                           float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // open, double, min
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val);
        }

        // open, float, min
        {
            adhoc::centered_canonical_open<RNGtest<std::uint64_t, min, max>,
                                           float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val);
        }

        // closed, double, max
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, float, max
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>,
                                             float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, double, min
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // closed, float, min
        {
            adhoc::centered_canonical_closed<RNGtest<std::uint64_t, min, max>,
                                             float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // default, double, max
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // default, float, max
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>,
                                              float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val);
        }

        // default, double, min
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // default, float, min
        {
            adhoc::centered_canonical_default<RNGtest<std::uint64_t, min, max>,
                                              float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }
    }

#ifndef _MSC_VER
    // 0, max, __uint128_t
    constexpr double small_val_128 = 2.9387358770557188e-39;
    {
        constexpr __uint128_t min = 0;
        constexpr __uint128_t max = std::numeric_limits<__uint128_t>::max();

        // open, double, max
        {
            adhoc::centered_canonical_open<RNGtest<__uint128_t, min, max>> gen(
                max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val_128);
        }

        // open, float, max
        {
            adhoc::centered_canonical_open<RNGtest<__uint128_t, min, max>,
                                           float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, -small_val_128);
        }

        // open, double, min
        {
            adhoc::centered_canonical_open<RNGtest<__uint128_t, min, max>> gen(
                min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val_128);
        }

        // open, float, min
        {
            adhoc::centered_canonical_open<RNGtest<__uint128_t, min, max>,
                                           float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, small_val_128);
        }

        // closed, double, max
        {
            adhoc::centered_canonical_closed<RNGtest<__uint128_t, min, max>>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, float, max
        {
            adhoc::centered_canonical_closed<RNGtest<__uint128_t, min, max>,
                                             float>
                gen(max);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), -1.0);
        }

        // closed, double, min
        {
            adhoc::centered_canonical_closed<RNGtest<__uint128_t, min, max>>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }

        // closed, float, min
        {
            adhoc::centered_canonical_closed<RNGtest<__uint128_t, min, max>,
                                             float>
                gen(min);

            auto val = gen();
            EXPECT_EQUAL(val, 0);
            EXPECT_EQUAL(std::copysign(1.0, val), 1.0);
        }
    }
#endif

    // LWG issue 2524 not there!
    {
        adhoc::exponential_distribution<float> dis(1.0);
        RNGtest<std::uint64_t, 0, std::numeric_limits<std::uint64_t>::max()>
            rng(std::numeric_limits<std::uint64_t>::max());

        auto val = dis(rng);
        EXPECT_NOT_EQUAL(val, std::numeric_limits<float>::infinity());
    }

    TEST_END;
}
