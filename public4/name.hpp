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

#ifndef ADHOC4_NAME_HPP
#define ADHOC4_NAME_HPP

#include "adhoc.hpp"
#include "differential_operator.hpp"
#include "utils/array_char.hpp"

#include <algorithm>
#include <array>
#include <string_view>

namespace adhoc4 {

template <detail::StringLiteral literal>
inline constexpr auto name(double_t<literal> /* t */) {
    constexpr auto tmp = std::to_array(literal.value);

    // we remove the last null character
    std::array<char, tmp.size() - 1> res;
    std::copy_n(tmp.data(), res.size(), res.data());
    return res;
}

#define AD_DEFINE_UNIOP_NAME(TYPE, OPNAME)                                     \
    template <class Input>                                                     \
    inline constexpr auto opname(Base<TYPE<Input>> /* t */) {                  \
        constexpr std::string_view tmp{#OPNAME};                               \
        std::array<char, tmp.size()> res;                                      \
        std::copy_n(tmp.data(), tmp.size(), res.data());                       \
        return res;                                                            \
    }

// Exponential functions
AD_DEFINE_UNIOP_NAME(exp_t, exp)
AD_DEFINE_UNIOP_NAME(log_t, log)

// Power functions
AD_DEFINE_UNIOP_NAME(sqrt_t, sqrt)

// Trigonometric functions
AD_DEFINE_UNIOP_NAME(sin_t, sin)
AD_DEFINE_UNIOP_NAME(cos_t, cos)
AD_DEFINE_UNIOP_NAME(tan_t, tan)
AD_DEFINE_UNIOP_NAME(asin_t, asin)
AD_DEFINE_UNIOP_NAME(acos_t, acos)
AD_DEFINE_UNIOP_NAME(atan_t, atan)

// Hyperbolic functions
AD_DEFINE_UNIOP_NAME(sinh_t, sinh)
AD_DEFINE_UNIOP_NAME(cosh_t, cosh)
AD_DEFINE_UNIOP_NAME(tanh_t, tanh)
AD_DEFINE_UNIOP_NAME(asinh_t, asinh)
AD_DEFINE_UNIOP_NAME(acosh_t, acosh)
AD_DEFINE_UNIOP_NAME(atanh_t, atanh)

// Error and gamma functions
AD_DEFINE_UNIOP_NAME(erfc_t, erfc)
AD_DEFINE_UNIOP_NAME(tgamma_t, tgamma)
AD_DEFINE_UNIOP_NAME(lgamma_t, lgamma)

// Mathematical special functions
AD_DEFINE_UNIOP_NAME(comp_ellint_1_t, comp_ellint_1)
AD_DEFINE_UNIOP_NAME(comp_ellint_2_t, comp_ellint_2)
AD_DEFINE_UNIOP_NAME(riemann_zeta_t, riemann_zeta)

template <template <class> class Univariate, class Input>
inline constexpr auto name(Base<Univariate<Input>> t) {
    return concatenate(opname(t), std::array<char, 1>{'('}, name(Input{}),
                       std::array<char, 1>{')'});
}

template <class Input1, class Input2>
inline constexpr auto opnamebi(Base<add_t<Input1, Input2>> /* t */) {
    constexpr std::string_view tmp{"+"};
    std::array<char, tmp.size()> res;
    std::copy_n(tmp.data(), tmp.size(), res.data());
    return res;
}

template <class Input1, class Input2>
inline constexpr auto opnamebi(Base<sub_t<Input1, Input2>> /* t */) {
    constexpr std::string_view tmp{"-"};
    std::array<char, tmp.size()> res;
    std::copy_n(tmp.data(), tmp.size(), res.data());
    return res;
}

template <class Input1, class Input2>
inline constexpr auto opnamebi(Base<mul_t<Input1, Input2>> /* t */) {
    constexpr std::string_view tmp{"*"};
    std::array<char, tmp.size()> res;
    std::copy_n(tmp.data(), tmp.size(), res.data());
    return res;
}

template <class Input1, class Input2>
inline constexpr auto opnamebi(Base<div_t<Input1, Input2>> /* t */) {
    constexpr std::string_view tmp{"/"};
    std::array<char, tmp.size()> res;
    std::copy_n(tmp.data(), tmp.size(), res.data());
    return res;
}

template <template <class, class> class Bivariate, class Input1, class Input2>
inline constexpr auto name(Base<Bivariate<Input1, Input2>> t) {
    return concatenate(std::array<char, 1>{'('}, name(Input1{}), opnamebi(t),
                       name(Input2{}), std::array<char, 1>{')'});
}

namespace detail {

inline constexpr auto digits(std::size_t val) {
    return static_cast<std::size_t>(std::trunc(std::log10(val))) + 1;
}

inline constexpr auto lastDigit(std::size_t val) { return (val % 10); }

inline constexpr auto single_digits(std::size_t val) {
    constexpr char digits[] = "0123456789";
    return digits[val];
}

template <std::size_t Val, std::size_t Idx = 0, std::size_t Size>
inline constexpr auto to_array_char_aux(std::array<char, Size> &result) {
    constexpr auto lastdigit = lastDigit(Val);
    result[Size - 1 - Idx] = single_digits(lastdigit);
    if constexpr (Idx < (Size - 1)) {
        constexpr auto remainder = (Val - lastdigit) / 10;
        to_array_char_aux<remainder, Idx + 1>(result);
    }
}

template <std::size_t Val> inline constexpr auto to_array_char() {
    constexpr auto Size = digits(Val);
    std::array<char, Size> result;
    to_array_char_aux<Val>(result);
    return result;
}

template <std::size_t Order, class Id>
inline constexpr auto name_aux(der::d<Order, Id> /* var */) {
    constexpr auto in_var = name(Id{});
    if constexpr (in_var.front() == '(' && in_var.back() == ')') {
        if constexpr (Order == 1) {
            return concatenate(std::array<char, 1>{'d'}, name(Id{}));

        } else {
            return concatenate(std::array<char, 2>{'d', '^'},
                               detail::to_array_char<Order>(), name(Id{}));
        }
    } else {
        if constexpr (Order == 1) {
            return concatenate(std::array<char, 1>{'d'},
                               std::array<char, 1>{'('}, name(Id{}),
                               std::array<char, 1>{')'});

        } else {
            return concatenate(
                std::array<char, 2>{'d', '^'}, detail::to_array_char<Order>(),
                std::array<char, 1>{'('}, name(Id{}), std::array<char, 1>{')'});
        }
    }
}

} // namespace detail

template <std::size_t Order1, class Id1, std::size_t... Order, class... Id>
inline constexpr auto
name(std::tuple<der::d<Order1, Id1>, der::d<Order, Id>...> /* t */) {

    constexpr auto thisvar = detail::name_aux(der::d<Order1, Id1>{});

    if constexpr (sizeof...(Id)) {
        return concatenate(thisvar, std::array<char, 1>{'*'},
                           name(std::tuple<der::d<Order, Id>...>{}));

    } else {
        return thisvar;
    }
}

} // namespace adhoc4

#endif // ADHOC4_NAME_HPP
