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

#ifndef ADHOC4_UNIVARIATE_TAN_HPP
#define ADHOC4_UNIVARIATE_TAN_HPP

#include "../base.hpp"
#include "../combinatorics/pascal.hpp"
#include "../utils/index_sequence.hpp"

#include <array>
#include <cmath>

namespace adhoc4 {

namespace detail {

template <std::size_t Idx, class IndexSequence>
constexpr auto Multiply_aux(double &result, double tan2, IndexSequence i) {
    result *= tan2;
    result += static_cast<double>(get<Idx - 1>(i));

    if constexpr (Idx > 2) {
        Multiply_aux<Idx - 2>(result, tan2, i);
    }
}

template <int First, int... I>
inline auto Multiply(double tan, double tan2,
                     std::integer_sequence<int, First, I...> i) {
    auto result = static_cast<double>(get<i.size() - 1>(i));
    if constexpr (First != 0) {
        Multiply_aux<i.size() - 2>(result, tan2, i);
    } else {
        Multiply_aux<i.size() - 2>(result, tan2, i);
        result *= tan;
    }
    return result;
}

template <bool Diff, std::size_t Idx, std::size_t Order, std::size_t Output,
          class IndexSequence>
inline void tan_aux(double tan, double tan2, std::array<double, Output> &res,
                    IndexSequence i) {
    constexpr auto next_i = NextPascal<Diff>(i);
    res[Idx] = Multiply(tan, tan2, next_i);
    if constexpr (Idx < Order) {
        tan_aux<Diff, Idx + 1, Order>(tan, tan2, res, next_i);
    }
}

} // namespace detail

template <class Input> struct tan_t : public Base<tan_t<Input>> {
    static inline auto v(double in) -> double { return std::tan(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double /* in */,
                         std::array<double, Output> &res) {
        // we use f'(x) - (f(x))^2 - 1 = 0
        static_assert(Order <= Output);

        const double tan2 = thisv * thisv;

        // equivalent to
        // detail::tan_aux<false, 0, Order>(thisv, tan2, res,
        //                                  std::integer_sequence<int, 0, 1>{});

        if constexpr (Order >= 1) {
            res[0] = 1. + tan2;
        }

        if constexpr (Order >= 2) {
            detail::tan_aux<false, 1, Order>(
                thisv, tan2, res, std::integer_sequence<int, 1, 0, 1>{});
        }
    }
};

template <class Derived> auto tan(Base<Derived> /* in */) {
    return tan_t<Derived>{};
}

template <class Input> struct tanh_t : public Base<tanh_t<Input>> {
    static inline auto v(double in) -> double { return std::tanh(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double /* in */,
                         std::array<double, Output> &res) {
        // we use f'(x) + (f(x))^2 - 1 = 0
        static_assert(Order <= Output);

        double const tan2 = thisv * thisv;

        // equivalent to
        // detail::tan_aux<true, 0, Order>(thisv, tan2, res,
        //                                 std::integer_sequence<int, 0, 1>{});

        if constexpr (Order >= 1) {
            res[0] = 1. - tan2;
        }

        if constexpr (Order >= 2) {
            detail::tan_aux<true, 1, Order>(
                thisv, tan2, res, std::integer_sequence<int, 1, 0, -1>{});
        }
    }
};

template <class Derived> auto tanh(Base<Derived> /* in */) {
    return tanh_t<Derived>{};
}

} // namespace adhoc4

#endif // ADHOC4_UNIVARIATE_TAN_HPP
