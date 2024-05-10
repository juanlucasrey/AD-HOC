/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ADHOC4_BASE_HPP
#define ADHOC4_BASE_HPP

#include <array>
#include <cstddef>

namespace adhoc4 {

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

    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double /* in */,
                         std::array<double, Output> &res) {
        static_assert(Order <= Output);

        if constexpr (Order >= 1) {
            res[0] = -thisv * thisv;
        }

        if constexpr (Order >= 2) {
            detail::invders<1>(res, thisv);
        }
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

} // namespace adhoc4

#endif // ADHOC4_BASE_HPP