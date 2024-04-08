#ifndef ADHOC3_BASE_HPP
#define ADHOC3_BASE_HPP

#include <array>
#include <cstddef>

namespace adhoc3 {

template <class Input1, class Input2> struct add_t;
template <class Input1, class Input2> struct sub_t;
template <class Input1, class Input2> struct mul_t;
template <class Input1, class Input2> struct div_t;

template <class Derived> struct Base;

namespace detail {

template <std::size_t N, std::size_t Order>
inline void invders(std::array<double, Order> &res, double val) {
    res[N] = -static_cast<double>(N + 1) * res[N - 1] * val;
    if constexpr ((N + 1) < Order) {
        invders<N + 1>(res, val);
    }
}

} // namespace detail

template <class Input> struct inv_t : public Base<inv_t<Input>> {
    static inline auto v(double in) -> double { return 1.0 / in; }
    static inline auto d(double thisv, double /* in */) -> double {
        return -thisv * thisv;
    }

    template <std::size_t Order>
    static inline auto d2(double thisv, double /* in */)
        -> std::array<double, Order> {
        std::array<double, Order> res;

        if constexpr (Order >= 1) {
            // detail::invders<1>(res, thisv);
            res[0] = -thisv * thisv;
        }

        if constexpr (Order >= 2) {
            detail::invders<1>(res, thisv);
        }

        return res;
    }
};

// it would be nice to add operator to constexpr float types in this class
// however this proposal would have to pass
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p1045r1.html
// https://github.com/cplusplus/papers/issues/603
template <class Derived> struct Base {
    template <class Derived2> auto operator+(Base<Derived2> /* rhs */) const {
        return add_t<Derived, Derived2>{};
    }

    template <class Derived2> auto operator-(Base<Derived2> /* rhs */) const {
        return sub_t<Derived, Derived2>{};
    }

    template <class Derived2> auto operator*(Base<Derived2> /* rhs */) const {
        return mul_t<Derived, Derived2>{};
    }

    // template <class Derived2> auto operator/(Base<Derived2> /* rhs */) const
    // {
    //     return div_t<Derived, Derived2>{};
    // }
    template <class Derived2> auto operator/(Base<Derived2> /* rhs */) const {
        return mul_t<Derived, inv_t<Derived2>>{};
    }
};

} // namespace adhoc3

#endif // ADHOC_BASE3_HPP
