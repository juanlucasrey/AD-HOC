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

#ifndef DEPENDENCY_HPP
#define DEPENDENCY_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"

namespace adhoc4 {

template <class I1, class I2>
constexpr auto operator==(I1 /* lhs */, I2 /* rhs */) -> bool {
    return std::is_same_v<I1, I2>;
}

// forward declaration so that mult_t operator== calls add_t operator==
template <class I1, class I2, class I3, class I4>
constexpr auto operator==(add_t<I1, I2> /* lhs */, add_t<I3, I4> /* rhs */)
    -> bool;

template <class I1, class I2, class I3, class I4>
constexpr auto operator==(mul_t<I1, I2> /* lhs */, mul_t<I3, I4> /* rhs */)
    -> bool {
    // multiplication is symetrical
    return ((I1{} == I3{}) && (I2{} == I4{})) ||
           ((I1{} == I4{}) && (I2{} == I3{}));
}

template <class I1, class I2, class I3, class I4>
constexpr auto operator==(add_t<I1, I2> /* lhs */, add_t<I3, I4> /* rhs */)
    -> bool {
    // addition is symetrical
    return ((I1{} == I3{}) && (I2{} == I4{})) ||
           ((I1{} == I4{}) && (I2{} == I3{}));
}

template <constants::ArgType N, class Denominator>
constexpr auto depends(constants::CD<N> /* num */, Denominator /* den */)
    -> bool {
    // why false? Because a constant does not depend on any variable
    return false;
};

template <detail::StringLiteral literal1, class Denominator>
constexpr auto depends(double_t<literal1> num, Denominator den) -> bool {
    return num == den;
};

template <template <class...> class Xvariate, class... Inputs,
          class Denominator>
constexpr auto depends(Xvariate<Inputs...> num, Denominator den) -> bool {
    if constexpr (num == den) {
        return true;
    } else {
        return ((depends(Inputs{}, den)) || ...);
    }
};

template <class I> constexpr auto is_input(I /* in */) -> bool { return false; }

template <detail::StringLiteral literal>
constexpr auto is_input(double_t<literal> /* in */) -> bool {
    return true;
}

template <class I> constexpr auto is_constant(I /* in */) -> bool {
    return false;
}

template <constants::ArgType N>
constexpr auto is_constant(constants::CD<N> /* in */) -> bool {
    return true;
}

} // namespace adhoc4

#endif // DEPENDENCY_HPP
