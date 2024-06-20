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

#include "../tests4/type_name.hpp"

#include <algorithm>
#include <array>

namespace adhoc4 {

template <detail::StringLiteral literal>
inline constexpr auto name(double_t<literal> /* t */) {
    constexpr auto tmp = std::to_array(literal.value);
    std::array<char, tmp.size() - 1> res;
    std::copy_n(tmp.data(), res.size(), res.data());
    return res;
}

#define AD_DEFINE_UNIOP_NAME(TYPE, OPNAME)                                     \
    template <class Input> inline constexpr auto opname(TYPE<Input> /* t */) { \
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
inline constexpr auto name(Univariate<Input> t) {
    constexpr auto zero = opname(t);
    constexpr std::array<char, 1> one{'('};
    constexpr auto two = name(Input{});
    constexpr std::array<char, 1> three{')'};
    constexpr unsigned N = zero.size() + one.size() + two.size() + three.size();
    std::array<char, N> result = {};

    auto it = result.begin();
    it = std::copy_n(zero.data(), zero.size(), it);
    it = std::copy_n(one.data(), one.size(), it);
    it = std::copy_n(two.data(), two.size(), it);
    it = std::copy_n(three.data(), three.size(), it);
    return result;
}

} // namespace adhoc4

#endif // ADHOC4_NAME_HPP
