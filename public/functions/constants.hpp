#ifndef FUNCTIONS_CONSTANTS_HPP
#define FUNCTIONS_CONSTANTS_HPP

#include <limits>

namespace adhoc2::constants {

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

} // namespace detail

template <class T> auto constexpr pi() -> T {
    return detail::bailey_borwein_plouffe<T>(0., 1., 0., 1.);
}

template <class T> auto constexpr sqrt(T x) -> T {
    return x >= 0 && x < std::numeric_limits<T>::infinity()
               ? detail::sqrtNewtonRaphson<T>(x, x, 0)
               : std::numeric_limits<T>::quiet_NaN();
}

} // namespace adhoc2::constants

#endif // FUNCTIONS_CONSTANTS_HPP
