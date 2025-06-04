#include "../../test_simple/test_simple_include.hpp"

#include <distribution/centered_canonical_distribution.hpp>

#include <cstdint>
#include <iostream>

template <class UIntType, UIntType Min, UIntType Max> struct RNGtest final {

    using result_type = UIntType;

    explicit RNGtest(result_type value_in) : value(value_in) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        return value;
    }

    static constexpr auto min() -> UIntType { return Min; }
    static constexpr auto max() -> UIntType { return Max; }

  private:
    UIntType value = Min;
};

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

    TEST_END;
}
