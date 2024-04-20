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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_DERIVATIVE_NON_NULL_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_DERIVATIVE_NON_NULL_HPP

#include "../adhoc.hpp"
#include "differential_operator.hpp"

namespace adhoc3 {

namespace detail {

template <std::size_t... Args> auto constexpr sum() -> std::size_t {
    return (Args + ...);
}

template <std::size_t Input1, std::size_t Input2, std::size_t... Inputs>
auto constexpr max() -> std::size_t {
    auto constexpr thismax = Input1 > Input2 ? Input1 : Input2;
    if constexpr (sizeof...(Inputs) == 0) {
        return thismax;
    } else {
        return max<thismax, Inputs...>();
    }
}

template <std::size_t Input1> auto constexpr max() -> std::size_t {
    return Input1;
}

// once we hit max, which represents infinity, we stay there
template <std::size_t Input1, std::size_t... Inputs>
auto constexpr sum_no_overflow2() -> std::size_t {
    if constexpr (sizeof...(Inputs) == 0) {
        return Input1;
    } else {
        if (std::numeric_limits<std::size_t>::max() - Input1 <
            sum_no_overflow<Inputs...>()) {
            return std::numeric_limits<std::size_t>::max();
        }

        if (std::numeric_limits<std::size_t>::max() -
                sum_no_overflow<Inputs...>() <
            Input1) {
            return std::numeric_limits<std::size_t>::max();
        }

        return Input1 + sum_no_overflow<Inputs...>();
    }
}

template <class IdNum, std::size_t OrderNum, std::size_t PowerNum,
          std::size_t N, std::size_t PowerDen>
auto constexpr derivative_non_null_simple(
    std::tuple<der2::p<PowerNum, der2::d<OrderNum, IdNum>>> /* ids1 */,
    der2::p<PowerDen, der2::d<1, double_t<N>>> /* id2 */) -> bool {
    auto constexpr summation = order<IdNum, double_t<N>>::call() * PowerNum;
    return summation >= PowerDen;
}

} // namespace detail

template <class... Ids, std::size_t... Orders1, std::size_t... Powers1,
          std::size_t... Ns, std::size_t... Orders2, std::size_t... Powers2>
// constexpr auto derivative_non_null(
auto derivative_non_null(
    std::tuple<der2::p<Powers1, der2::d<Orders1, Ids>>...> /* ids1 */,
    std::tuple<
        der2::p<Powers2, der2::d<Orders2, double_t<Ns>>>...> /* ids2 */) {

    auto constexpr numer_pow = detail::sum<Orders1 * Powers1...>();
    auto constexpr denom_pow = detail::sum<Orders2 * Powers2...>();

    if constexpr (denom_pow != numer_pow) {
        return false;
    } else {
        auto constexpr numer_max_d = detail::max<Orders1...>();
        auto constexpr denom_max_d = detail::max<Orders2...>();
        if constexpr (denom_max_d > numer_max_d) {
            return false;
        } else {
            if (denom_max_d == 1 && sizeof...(Ids) == 1) {
                return (
                    detail::derivative_non_null_simple(
                        std::tuple<
                            der2::p<Powers1, der2::d<Orders1, Ids>>...>{},
                        der2::p<Powers2, der2::d<Orders2, double_t<Ns>>>{}) &&
                    ...);
            } else {
                return true;
            }
        }
    }
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_DERIVATIVE_NON_NULL_HPP
