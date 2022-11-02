#ifndef CONSTANTS_DOUBLE_HPP
#define CONSTANTS_DOUBLE_HPP

#include <base.hpp>
#include <constants_constexpr.hpp>
#include <constants_type.hpp>

#include <array>
#include <cstdint>

namespace adhoc2::constants {

namespace detail {

template <bool positive, unsigned int End, unsigned int Count = 1>
auto constexpr uint64_to_double_aux2(double x, std::uint32_t e) -> double {

    constexpr std::array<double, 10> pow2{
        2.,                     // 2^0b1
        4.,                     // 2^0b10
        16.,                    // 2^0b100
        256.,                   // 2^0b1000
        65536.,                 // 2^0b10000
        4294967296.,            // 2^0b100000
        1.8446744073709552e+19, // 2^0b1000000
        3.4028236692093846e+38, // 2^0b10000000
        1.157920892373162e+77,  // 2^0b10000000
        1.3407807929942597e+154 // 2^0b100000000
    };

    if constexpr (Count == End) {
        return x;
    } else {
        bool apply = (e >> (Count - 1)) & 1U;
        if constexpr (positive) {
            return uint64_to_double_aux2<positive, End, Count + 1>(
                apply ? x * pow2[Count - 1] : x, e);
        } else {
            return uint64_to_double_aux2<positive, End, Count + 1>(
                apply ? x / pow2[Count - 1] : x, e);
        }
    }
}

auto constexpr uint64_to_double_aux(double x, std::int32_t e) -> double {
    // std::abs breaks constexpr!!
    // auto eu = static_cast<std::uint32_t>(std::abs(e));
    std::uint32_t eu =
        e < 0 ? static_cast<std::uint32_t>(-e) : static_cast<std::uint32_t>(e);

    if (e == 0U) {
        return x;
    } else if (e < 0) {
        // 11 is the size of the exponent. with this configuration we
        // are guaranteed 11 recursions only
        return uint64_to_double_aux2<false, 11U>(x, eu);
    } else {
        return uint64_to_double_aux2<true, 11U>(x, eu);
    }
}

auto constexpr uint64_to_double(std::uint64_t x) -> double {
    // based on
    // https://en.wikipedia.org/wiki/Double-precision_floating-point_format#Exponent_encoding

    bool sign = static_cast<bool>(x & 0x8000000000000000);
    std::uint64_t exp = (x >> 52U) & 0x7ffU;
    std::uint64_t val = (x & 0xfffffffffffffU);

    // exceptions inf and nan
    if (exp == 2047) {
        if (val == 0) {
            return sign ? -std::numeric_limits<double>::infinity()
                        : std::numeric_limits<double>::infinity();
        }
        return sign ? -std::numeric_limits<double>::quiet_NaN()
                    : std::numeric_limits<double>::quiet_NaN();
    }

    // subnormal case
    if (exp == 0) {
        if (val == 0) {
            return sign ? -0 : 0;
        }
        auto vald = static_cast<double>(val);
        double unsigned_result = vald * 4.9406564584124654e-324;
        return sign ? -unsigned_result : unsigned_result;
    }

    val |= 0x10000000000000U;

    auto exp_trunc = static_cast<std::int32_t>(exp) - 0x3ff;
    auto denom = static_cast<double>(0x10000000000000U);
    auto vald = static_cast<double>(val);

    double unsigned_result = uint64_to_double_aux(vald / denom, exp_trunc);
    return sign ? -unsigned_result : unsigned_result;
}

constexpr std::uint64_t double_to_uint64(double half) {
    if (half > 0.5) {
        return 1;
    } else {
        return 2;
    }
}

} // namespace detail

template <std::uint64_t D>
struct Double : public Base<Double<D>>, public ConstBase<Double<D>> {
    constexpr auto static v() -> double { return detail::uint64_to_double(D); }
};

} // namespace adhoc2::constants

#endif // CONSTANTS_DOUBLE_HPP
