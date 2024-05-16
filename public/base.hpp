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

#ifndef ADHOC_BASE_HPP
#define ADHOC_BASE_HPP

#include <cstddef>

namespace adhoc {

template <class Input1, class Input2> struct add_t;
template <class Input1, class Input2> struct sub_t;
template <class Input1, class Input2> struct mul_t;
template <class Input1, class Input2> struct div_t;

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

    template <class Derived2> auto operator/(Base<Derived2> /* rhs */) const {
        return div_t<Derived, Derived2>{};
    }
};

} // namespace adhoc

#endif // ADHOC_BASE_HPP
