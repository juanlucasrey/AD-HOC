#ifndef CONSTANTS_DOUBLE_HPP
#define CONSTANTS_DOUBLE_HPP

#include <base.hpp>
#include <constants_constexpr.hpp>
#include <constants_type.hpp>

namespace adhoc2::constants {

namespace detail {

template <bool positive, unsigned int End, unsigned int Count = 1>
auto constexpr uint64_to_double_aux2(double x, std::uint32_t e,
                                     double multiplier = 2.0) -> double {
    if constexpr (Count == End) {
        return x;
    } else {
        if constexpr (positive) {
            bool apply = (e >> (Count - 1)) & 1U;
            return uint64_to_double_aux2<positive, End, Count + 1>(
                apply ? x * multiplier : x, e, multiplier * multiplier);
        } else {
            bool apply = (e >> (Count - 1)) & 1U;
            return uint64_to_double_aux2<positive, End, Count + 1>(
                apply ? x / multiplier : x, e, multiplier * multiplier);
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
