#ifndef CONSTANTS_TYPE_HPP
#define CONSTANTS_TYPE_HPP

#include <base.hpp>
#include <constants_constexpr.hpp>

namespace adhoc2::constants {

template <class Derived> struct ConstBase {
    operator double() const { return Derived::v(); }
};

struct Pi : public Base<Pi>, public ConstBase<Pi> {
    constexpr auto static v() -> double { return detail::pi<double>(); }
};

template <std::size_t N>
struct Const : public Base<Const<N>>, public ConstBase<Const<N>> {
    constexpr auto static v() -> double { return static_cast<double>(N); }
};

template <class N, class D>
struct Frac : public Base<Frac<N, D>>, public ConstBase<Frac<N, D>> {
    constexpr auto static v() -> double { return N::v() / D::v(); }
};

template <class V>
struct Sqrt : public Base<Sqrt<V>>, public ConstBase<Sqrt<V>> {
    constexpr auto static v() -> double { return detail::sqrt(V::v()); }
};

template <class V>
struct Minus : public Base<Minus<V>>, public ConstBase<Minus<V>> {
    constexpr auto static v() -> double { return -V::v(); }
};

} // namespace adhoc2::constants

#endif // CONSTANTS_TYPE_HPP
