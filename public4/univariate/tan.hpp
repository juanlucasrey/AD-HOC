/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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

template <bool PrevPlus, bool NextPlus, std::size_t Order, std::size_t N = 0,
          std::size_t Output, class IndexSequence>
inline void tan_aux(double tan, double tan2, std::array<double, Output> &res,
                    IndexSequence i) {
    constexpr auto next_i = NextPascal<PrevPlus, NextPlus>(i);
    res[N] = Multiply(tan, tan2, next_i);
    if constexpr ((N + 1) < Order) {
        tan_aux<PrevPlus, NextPlus, Order, N + 1>(tan, tan2, res, next_i);
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
        detail::tan_aux<true, true, Order>(thisv, tan2, res,
                                           std::integer_sequence<int, 0, 1>{});
    }
};

template <class Derived> constexpr auto tan(Base<Derived> /* in */) {
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
        detail::tan_aux<false, true, Order>(thisv, tan2, res,
                                            std::integer_sequence<int, 0, 1>{});
    }
};

template <class Derived> constexpr auto tanh(Base<Derived> /* in */) {
    return tanh_t<Derived>{};
}

} // namespace adhoc4

#endif // ADHOC4_UNIVARIATE_TAN_HPP
