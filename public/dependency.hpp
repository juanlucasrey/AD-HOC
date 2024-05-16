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

#ifndef ADHOC_DEPENDENCY_HPP
#define ADHOC_DEPENDENCY_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"

#include <type_traits>

namespace adhoc {

template <class First, class Second> struct equal_or_depends;

template <constants::ArgType N, class Second>
struct equal_or_depends<constants::CD<N>, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<constants::CD<N>, Second>;
    }
};

template <std::size_t N, class Second>
struct equal_or_depends<double_t<N>, Second> {
    constexpr static auto call() noexcept -> bool {
        return std::is_same_v<double_t<N>, Second>;
    }
};

template <template <class> class Univariate, class Input, class Second>
struct equal_or_depends<Univariate<Input>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(
            std::is_convertible_v<Univariate<Input>, Base<Univariate<Input>>>);
        return std::is_same_v<Univariate<Input>, Second> ||
               equal_or_depends<Input, Second>::call();
    }
};

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
struct equal_or_depends<Bivariate<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(std::is_convertible_v<Bivariate<Input1, Input2>,
                                            Base<Bivariate<Input1, Input2>>>);
        return std::is_same_v<Bivariate<Input1, Input2>, Second> ||
               equal_or_depends<Input1, Second>::call() ||
               equal_or_depends<Input2, Second>::call();
    }
};

template <class First, class Second> struct depends;

template <constants::ArgType N, class Second>
struct depends<constants::CD<N>, Second> {
    constexpr static auto call() noexcept -> bool { return false; }
};

template <std::size_t N, class Second> struct depends<double_t<N>, Second> {
    constexpr static auto call() noexcept -> bool { return false; }
};

template <template <class> class Univariate, class Input, class Second>
struct depends<Univariate<Input>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(
            std::is_convertible_v<Univariate<Input>, Base<Univariate<Input>>>);
        return equal_or_depends<Input, Second>::call();
    }
};

template <template <class, class> class Bivariate, class Input1, class Input2,
          class Second>
struct depends<Bivariate<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> bool {
        static_assert(std::is_convertible_v<Bivariate<Input1, Input2>,
                                            Base<Bivariate<Input1, Input2>>>);
        return equal_or_depends<Input1, Second>::call() ||
               equal_or_depends<Input2, Second>::call();
    }
};

template <class First, class... Others>
inline constexpr auto equal_or_depends_many() noexcept -> bool {
    if constexpr (sizeof...(Others) == 0) {
        return false;
    } else {
        return (equal_or_depends<First, Others>::call() || ...);
    }
}

template <class First, class... Others>
inline constexpr auto depends_many() noexcept -> bool {
    if constexpr (sizeof...(Others) == 0) {
        return false;
    } else {
        return (depends<First, Others>::call() || ...);
    }
}

} // namespace adhoc

#endif // ADHOC_DEPENDENCY_HPP
