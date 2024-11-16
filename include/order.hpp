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

#ifndef ORDER_HPP
#define ORDER_HPP

#include "add_sat.hpp"
#include "adhoc.hpp"
#include "constants_type.hpp"

namespace adhoc {

//  order(Numerator, Denominator) gives the highest order for which
//  d^n(First)/d(Second)^n is non-null
template <constants::ArgType N, class Denominator>
constexpr auto order(constants::CD<N> /* num */, Denominator /* den */)
    -> std::size_t {
    // why 0? Because Second has to be a variable and therefore could never be a
    // constant
    return 0;
};

template <detail::StringLiteral literal, class Denominator>
constexpr auto order(double_t<literal> /* num */, Denominator /* den */)
    -> std::size_t {
    // we can't drill down further here, so it's either the variable or it isn't
    return std::is_same_v<double_t<literal>, Denominator> ? 1UL : 0UL;
};

template <template <class> class Univariate, class Input, class Denominator>
constexpr auto order(Univariate<Input> /* num */, Denominator den)
    -> std::size_t {
    return std::is_same_v<Univariate<Input>, Denominator> ? 1UL
           : order(Input{}, den) ? std::numeric_limits<std::size_t>::max()
                                 : 0UL;
};

template <class Input1, class Input2, class Denominator>
constexpr auto order(add_t<Input1, Input2> /* num */, Denominator den)
    -> std::size_t {
    static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                        Base<add_t<Input1, Input2>>>);
    return std::is_same_v<add_t<Input1, Input2>, Denominator>
               ? 1UL
               : std::max(order(Input1{}, den), order(Input2{}, den));
};

template <class Input1, class Input2, class Denominator>
constexpr auto order(sub_t<Input1, Input2> /* num */, Denominator den)
    -> std::size_t {
    static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                        Base<add_t<Input1, Input2>>>);
    return std::is_same_v<add_t<Input1, Input2>, Denominator>
               ? 1UL
               : std::max(order(Input1{}, den), order(Input2{}, den));
};

template <class Input1, class Input2, class Denominator>
constexpr auto order(mul_t<Input1, Input2> /* num */, Denominator den)
    -> std::size_t {
    static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                        Base<add_t<Input1, Input2>>>);
    return std::is_same_v<add_t<Input1, Input2>, Denominator>
               ? 1UL
               : detail::add_sat(order(Input1{}, den), order(Input2{}, den));
};

template <class Input1, class Input2, class Denominator>
constexpr auto order(div_t<Input1, Input2> /* num */, Denominator den)
    -> std::size_t {
    static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                        Base<add_t<Input1, Input2>>>);

    // we treat a/b as a * inv(b), inv being a typical univariate operator.
    // since inv_t does not exist, we use exp_t which would return the same
    // order.
    return std::is_same_v<add_t<Input1, Input2>, Denominator>
               ? 1UL
               : detail::add_sat(order(Input1{}, den),
                                 order(exp_t<Input2>{}, den));
};

} // namespace adhoc

#endif // ORDER_HPP
