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

#ifndef UTILS_FLOAT_TO_STRING_HPP
#define UTILS_FLOAT_TO_STRING_HPP

#include "../constants_type.hpp"
#include "../name.hpp"
#include "array_char.hpp"

#include <array>
#include <cmath>

namespace adhoc {

namespace detail {
inline constexpr auto exponent(double val) {
    return static_cast<int>(std::floor(std::log10(val)));
}

template <int Val> inline constexpr auto to_array_char_int() {
    if constexpr (Val >= 0) {
        return to_array_char<Val>();
    } else {
        return concatenate(std::array<char, 1>{'-'}, to_array_char<-Val>());
    }
}

template <class C> inline constexpr auto exponent_string_aux(C val) {
    constexpr auto valin = val.v();
    constexpr auto expon = exponent(valin);
    if constexpr (expon == 0) {
        return std::array<char, 0>{};
    } else if constexpr (expon > 0) {
        constexpr auto arrayexp = to_array_char_int<expon>();
        return concatenate(std::array<char, 2>{'e', '+'}, arrayexp);
    } else {
        constexpr auto arrayexp = to_array_char_int<expon>();
        return concatenate(std::array<char, 1>{'e'}, arrayexp);
    }
}

inline constexpr auto first_digit(double val) {
    return static_cast<std::size_t>(
        std::floor(val * std::pow(10, -exponent(val))));
}

template <std::size_t End, std::size_t Pos, int PrevExp, class C1, class C2,
          class AS>
constexpr inline auto value_string_aux(C1 val, C2 mindiff, AS &result) {
    constexpr auto valin = val.v();
    constexpr auto expon_d = exponent(valin);
    if constexpr (expon_d == PrevExp - 1) {
        if constexpr ((Pos + 1) == End) {
            constexpr auto valin2 = valin * std::pow(10, -expon_d);
            constexpr std::size_t dig = std::round(valin2);
            // if constexpr (dig != 0) {
            result[Pos] = single_digits(dig);
            // }
        } else {
            constexpr auto first_d = first_digit(valin);
            result[Pos] = single_digits(first_d);
            if constexpr ((Pos + 1) < End) {
                constexpr auto newval = valin - static_cast<double>(first_d) *
                                                    std::pow(10, expon_d);

                constexpr auto mindiffin = mindiff.v();
                if constexpr (newval > mindiffin) {
                    value_string_aux<End, Pos + 1, PrevExp - 1>(
                        constants::CD<newval>(), mindiff, result);
                }
            }
        }

    } else {
        result[Pos] = '0';
        if constexpr ((Pos + 1) < End) {
            value_string_aux<End, Pos + 1, PrevExp - 1>(val, mindiff, result);
        }
    }
}

template <std::size_t MaxSize, std::size_t CurrentSize, int PrevExp, class C1,
          class C2>
constexpr inline auto value_string_aux_size2(C1 val, C2 mindiff) {
    constexpr auto valin = val.v();
    constexpr auto expon_d = exponent(valin);
    if constexpr (expon_d == PrevExp - 1) {
        constexpr auto first_d = first_digit(valin);
        if constexpr ((CurrentSize + 1) < MaxSize) {
            if constexpr ((CurrentSize + 1) == MaxSize) {
                constexpr auto valin2 = valin * std::pow(10, -expon_d);
                constexpr std::size_t dig = std::round(valin2);
                // if constexpr (dig != 0) {
                return CurrentSize + 1;
                // }
            } else {
                constexpr auto newval = valin - static_cast<double>(first_d) *
                                                    std::pow(10, expon_d);
                constexpr auto mindiffin = mindiff.v();
                if constexpr (newval > mindiffin) {
                    return value_string_aux_size2<MaxSize, CurrentSize + 1,
                                                  PrevExp - 1>(
                        constants::CD<newval>(), mindiff);
                } else {
                    return CurrentSize + 1;
                }
            }
        } else {
            return CurrentSize;
        }
    } else {
        if constexpr ((CurrentSize + 1) < MaxSize) {
            return value_string_aux_size2<MaxSize, CurrentSize + 1,
                                          PrevExp - 1>(val, mindiff);
        } else {
            return CurrentSize + 1;
        }
    }
}

template <std::size_t MaxSize, class C1, class C2>
constexpr inline auto value_string_aux_size(C1 val, C2 mindiff) {
    constexpr auto valin = val.v();
    constexpr auto mindiffin = mindiff.v();
    constexpr auto first_d = first_digit(valin);
    constexpr auto expon_d = exponent(valin);
    constexpr auto valin2 =
        valin - static_cast<double>(first_d) * std::pow(10, expon_d);
    if constexpr (valin2 > mindiffin) {
        return value_string_aux_size2<MaxSize, 2, expon_d>(
            constants::CD<valin2>(), mindiff);
    } else {
        return 1;
    }
}

template <class C> constexpr inline auto value_string_aux(C val) {
    constexpr auto valin = val.v();

    constexpr auto min_diff =
        std::nextafter(valin, std::numeric_limits<double>::max()) - valin;

    constexpr std::size_t max_size = 19;
    constexpr std::size_t size = value_string_aux_size<max_size>(
        constants::CD<valin>(), constants::CD<min_diff>());

    std::array<char, size> result;
    result.fill('0');

    constexpr auto first_d = first_digit(valin);
    constexpr auto expon_d = exponent(valin);
    result[0] = single_digits(first_d);
    constexpr auto valin2 =
        valin - static_cast<double>(first_d) * std::pow(10, expon_d);

    if constexpr (valin2 > min_diff) {
        result[1] = '.';
        value_string_aux<result.size(), 2, expon_d>(
            constants::CD<valin2>(), constants::CD<min_diff>(), result);
    }

    return result;
}

template <class C> constexpr inline auto name_positive(C val) {
    constexpr auto fraction_string = value_string_aux(val);
    constexpr auto exponent_string = exponent_string_aux(val);
    return concatenate(fraction_string, exponent_string);
}

} // namespace detail

template <constants::ArgType N>
constexpr inline auto name(constants::CD<N> val) {
    constexpr auto valin = val.v();
    if constexpr (valin == 0) {
        return std::array<char, 1>{'0'};
    } else if constexpr (valin >= 0) {
        return detail::name_positive(val);
    } else {
        constexpr auto val_inv = constants::CD<-valin>();
        return concatenate(std::array<char, 1>{'-'},
                           detail::name_positive(val_inv));
    }
}

} // namespace adhoc

#endif // UTILS_FLOAT_TO_STRING_HPP
