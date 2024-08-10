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

#ifndef ADHOC4_ORDER_HPP
#define ADHOC4_ORDER_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"

namespace adhoc4 {

//  order<First, Second>::call() gives the highest order for which
//  d^n(First)/d(Second)^n is non-null
template <class First, class Second> struct order;

template <constants::ArgType N, class Second>
struct order<constants::CD<N>, Second> {
    // why 0? Because Second has to be a variable and therefore could never be a
    // constant
    constexpr static auto call() noexcept -> std::size_t { return 0; }
};

template <detail::StringLiteral literal, class Second>
struct order<double_t<literal>, Second> {
    // we can't drill down further here, so it's either the variable or it isn't
    constexpr static auto call() noexcept -> std::size_t {
        return std::is_same_v<double_t<literal>, Second> ? 1UL : 0UL;
    }
};

template <template <class> class Univariate, class Input, class Second>
struct order<Univariate<Input>, Second> {
    constexpr static auto call() noexcept -> std::size_t {
        return std::is_same_v<Univariate<Input>, Second> ? 1UL
               : order<Input, Second>::call()
                   ? std::numeric_limits<std::size_t>::max()
                   : 0UL;
    }
};

template <class Input1, class Input2, class Second>
struct order<add_t<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> std::size_t {
        static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                            Base<add_t<Input1, Input2>>>);
        return std::is_same_v<add_t<Input1, Input2>, Second>
                   ? 1UL
                   : std::max(order<Input1, Second>::call(),
                              order<Input2, Second>::call());
    }
};

template <class Input1, class Input2, class Second>
struct order<sub_t<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> std::size_t {
        static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                            Base<add_t<Input1, Input2>>>);
        return std::is_same_v<add_t<Input1, Input2>, Second>
                   ? 1UL
                   : std::max(order<Input1, Second>::call(),
                              order<Input2, Second>::call());
    }
};

namespace detail {

template <typename T, typename... Types>
concept is_all_same = (... && std::is_same<T, Types>::value);

constexpr auto add_sat(std::size_t lhs) { return lhs; }

template <is_all_same<std::size_t>... Types>
constexpr auto add_sat(std::size_t lhs, Types... args) {
    auto const rhs = add_sat(args...);

    if (std::numeric_limits<std::size_t>::max() - lhs < rhs) {
        return std::numeric_limits<std::size_t>::max();
    }

    if (std::numeric_limits<std::size_t>::max() - rhs < lhs) {
        return std::numeric_limits<std::size_t>::max();
    }

    return lhs + rhs;
}

} // namespace detail

template <class Input1, class Input2, class Second>
struct order<mul_t<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> std::size_t {
        static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                            Base<add_t<Input1, Input2>>>);
        return std::is_same_v<add_t<Input1, Input2>, Second>
                   ? 1UL
                   : detail::add_sat(order<Input1, Second>::call(),
                                     order<Input2, Second>::call());
    }
};

template <class Input1, class Input2, class Second>
struct order<div_t<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> std::size_t {
        static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                            Base<add_t<Input1, Input2>>>);

        // we treat a/b as a * inv(b), inv being a typical univariate operator.
        // since inv_t does not exist, we use exp_t which would return the same
        // order.
        return std::is_same_v<add_t<Input1, Input2>, Second>
                   ? 1UL
                   : detail::add_sat(order<Input1, Second>::call(),
                                     order<exp_t<Input2>, Second>::call());
    }
};

} // namespace adhoc4

#endif // ADHOC4_ORDER_HPP
