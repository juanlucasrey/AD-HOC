#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <limits>

#include <base.hpp>

namespace adhoc2::constants2 {

namespace detail {

template <class T>
auto constexpr bailey_borwein_plouffe(T curr, T prev, T N, T power) -> T {
    return curr == prev ? curr
                        : bailey_borwein_plouffe(
                              (1. / power *
                               (4. / (8. * N + 1.) - 2. / (8. * N + 4.) -
                                1. / (8. * N + 5.) - 1. / (8. * N + 6.))) +
                                  curr,
                              curr, N + 1., power * 16.);
}

template <class T> auto constexpr sqrtNewtonRaphson(T x, T curr, T prev) -> T {
    return curr == prev ? curr
                        : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}

template <class T> auto constexpr pi() -> T {
    return detail::bailey_borwein_plouffe<T>(0., 1., 0., 1.);
}

template <class T> auto constexpr sqrt(T x) -> T {
    return x >= 0 && x < std::numeric_limits<T>::infinity()
               ? detail::sqrtNewtonRaphson<T>(x, x, 0)
               : std::numeric_limits<T>::quiet_NaN();
}

} // namespace detail

template <class Derived> struct ConstBase {
    operator double() const { return Derived::v(); }
    template <class Denom> constexpr static auto depends() noexcept -> bool {
        return false;
    }
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

} // namespace adhoc2::constants2

#endif // CONSTANTS_HPP
