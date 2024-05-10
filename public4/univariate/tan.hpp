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

#include <array>
#include <cmath>

namespace adhoc4 {

namespace detail {

template <std::size_t N>
constexpr auto NextPascal(std::array<std::size_t, N> const &prev) {
    std::array<std::size_t, N + 1> arr;
    arr.fill(0);

    arr[0] = prev[1];

    for (std::size_t i = 1; i < (prev.size() - 1); i++) {
        arr[i] = prev[i - 1] * (i - 1) + prev[i + 1] * (i + 1);
    }

    arr.back() = prev.back() * (prev.size() - 1);
    return arr;
}

template <std::size_t... I, std::size_t... I2>
constexpr auto NextPascal_aux(std::index_sequence<I...> /* i */,
                              std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    constexpr auto next = NextPascal(temp);
    return std::index_sequence<next[I2]...>{};
}

template <std::size_t... I>
constexpr auto NextPascal(std::index_sequence<I...> i) {
    return NextPascal_aux(i, std::make_index_sequence<sizeof...(I) + 1>{});
}

template <std::size_t Idx, std::size_t... I>
constexpr auto get(std::index_sequence<I...>) {
    return std::array<std::size_t, sizeof...(I)>{I...}[Idx];
}

template <std::size_t Idx, std::size_t... I>
constexpr auto Multiply_aux(double &result, double tan2,
                            std::index_sequence<I...> i) {
    result *= tan2;
    result += static_cast<double>(get<Idx - 1>(i));

    if constexpr (Idx > 2) {
        Multiply_aux<Idx - 2>(result, tan2, i);
    }
}

template <std::size_t First, std::size_t... I>
inline auto Multiply(double tan, double tan2,
                     std::index_sequence<First, I...> i) {
    auto result = static_cast<double>(get<i.size() - 1>(i));
    if constexpr (First != 0) {
        Multiply_aux<i.size() - 2>(result, tan2, i);
    } else {
        Multiply_aux<i.size() - 2>(result, tan2, i);
        result *= tan;
    }
    return result;
}

template <std::size_t Idx, std::size_t Order, std::size_t Output,
          class IndexSequence>
inline void tan_aux(double tan, double tan2, std::array<double, Output> &res,
                    IndexSequence i) {
    constexpr auto next_i = NextPascal(i);
    res[Idx] = Multiply(tan, tan2, next_i);
    if constexpr (Idx < Order) {
        tan_aux<Idx + 1, Order>(tan, tan2, res, next_i);
    }
}

} // namespace detail

template <class Input> struct tan_t : public Base<tan_t<Input>> {
    static inline auto v(double in) -> double { return std::tan(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use f''(x) + f(x) = 0
        static_assert(Order <= Output);

        double const tan2 = thisv * thisv;

        // equivalent to
        // detail::tan_aux<0, Order>(thisv, tan2, res, std::index_sequence<0U,
        // 1U>{});

        if constexpr (Order >= 1) {
            res[0] = tan2 + 1;
        }

        if constexpr (Order >= 2) {
            detail::tan_aux<1, Order>(thisv, tan2, res,
                                      std::index_sequence<1U, 0U, 1U>{});
        }
    }
};

template <class Derived> auto tan(Base<Derived> /* in */) {
    return tan_t<Derived>{};
}

} // namespace adhoc4

#endif // ADHOC4_UNIVARIATE_TAN_HPP
