/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/adhoc2).
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

#ifndef ADHOC3_DEPENDENCY_HPP
#define ADHOC3_DEPENDENCY_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"

#include <type_traits>

namespace adhoc3 {

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

//  order<First, Second>::call() gives the highest order for which
//  d^n(First)/d(Second)^n is non-null
template <class First, class Second> struct order;

template <constants::ArgType N, class Second>
struct order<constants::CD<N>, Second> {
    // why 0? Because Second has to be a variable and therefore could never be a
    // constant
    constexpr static auto call() noexcept -> std::size_t { return 0; }
};

template <std::size_t N, class Second> struct order<double_t<N>, Second> {
    // we can't drill down further here, so it's either the variable or it isn't
    constexpr static auto call() noexcept -> std::size_t {
        return std::is_same_v<double_t<N>, Second> ? 1UL : 0UL;
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

// once we hit max, which represents infinity, we stay there
template <std::size_t Input1, std::size_t... Inputs>
auto constexpr sum_no_overflow() -> std::size_t {
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

} // namespace  detail

template <class Input1, class Input2, class Second>
struct order<mul_t<Input1, Input2>, Second> {
    constexpr static auto call() noexcept -> std::size_t {
        static_assert(std::is_convertible_v<add_t<Input1, Input2>,
                                            Base<add_t<Input1, Input2>>>);
        return std::is_same_v<add_t<Input1, Input2>, Second>
                   ? 1UL
                   : detail::sum_no_overflow<order<Input1, Second>::call(),
                                             order<Input2, Second>::call()>();
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
                   : detail::sum_no_overflow<
                         order<Input1, Second>::call(),
                         order<exp_t<Input2>, Second>::call()>();
    }
};

} // namespace adhoc3

#endif // ADHOC3_DEPENDENCY_HPP
