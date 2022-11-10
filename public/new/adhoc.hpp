#ifndef NEW_ADHOC_HPP
#define NEW_ADHOC_HPP

// #include <constants_constexpr.hpp>
// #include <filter.hpp>
#include "base.hpp"
// #include <utils.hpp>

// #include <array>
#include <cmath>
// #include <iostream>
// #include <type_traits>

namespace adhoc {

template <class Input> struct exp_t : public Base<exp_t<Input>> {
    static inline auto v(double in) -> double { return std::exp(in); }
    static inline auto d(double thisv, double /* inputv */) -> double {
        return thisv;
    }
};

template <class Derived>
auto exp(Base<Derived> const &) -> exp_t<const Derived> {
    return exp_t<const Derived>{};
}

template <class Input> struct cos_t : public Base<cos_t<Input>> {
    static inline auto v(double in) -> double { return std::cos(in); }
};

template <class Derived>
auto cos(Base<Derived> const &) -> cos_t<const Derived> {
    return cos_t<const Derived>{};
}

template <class Input> struct sqrt_t : public Base<sqrt_t<Input>> {
    static inline auto v(double in) -> double { return std::sqrt(in); }
};

template <class Derived>
auto sqrt(Base<Derived> const &) -> sqrt_t<const Derived> {
    return sqrt_t<const Derived>{};
}

template <class Input> struct log_t : public Base<log_t<Input>> {
    static inline auto v(double in) -> double { return std::log(in); }
};

template <class Derived>
auto log(Base<Derived> const &) -> log_t<const Derived> {
    return log_t<const Derived>{};
}

template <class Input> struct erfc_t : public Base<erfc_t<Input>> {
    static inline auto v(double in) -> double { return std::erfc(in); }
};

template <class Derived>
auto erfc(Base<Derived> const &) -> erfc_t<const Derived> {
    return erfc_t<const Derived>{};
}

template <class Input1, class Input2>
struct add_t : public Base<add_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 + in2; }
};

template <class Input1, class Input2>
struct sub_t : public Base<sub_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 - in2; }
};

template <class Input1, class Input2>
struct mul_t : public Base<mul_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 * in2; }
};

template <class Input1, class Input2>
struct div_t : public Base<div_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 / in2; }
};

namespace detail {
template <std::size_t N> struct adouble_aux : public Base<adouble_aux<N>> {};

} // namespace detail

} // namespace adhoc

#endif // NEW_ADHOC_HPP