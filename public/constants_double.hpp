#ifndef CONSTANTS_DOUBLE_HPP
#define CONSTANTS_DOUBLE_HPP

#include <base.hpp>
#include <constants_constexpr.hpp>
#include <constants_type.hpp>

namespace adhoc2::constants {

template <typename D>
struct Double : public Base<Double<D>>, public ConstBase<Double<D>> {
    constexpr auto static v() -> double { return 0.5; }
};

constexpr std::size_t convert_to_type(double half) {
    if (half > 0.5) {
        return 1;
    } else {
        return 2;
    }
    // constexpr std::size_t Numerator = 1;
    // constexpr std::size_t Denominator = 2;
    // return Frac<Const<Numerator>, Const<Denominator>>();
}

} // namespace adhoc2::constants

#endif // CONSTANTS_DOUBLE_HPP
