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

#ifndef BRNG_TOOLS_UNSHIFT
#define BRNG_TOOLS_UNSHIFT

#include "mask.hpp"

#include <algorithm>
#include <cstddef>

namespace adhoc {

namespace detail {

template <std::size_t End, std::size_t Shift, std::size_t Position = 0,
          class UIntType>
inline auto unshift_left_xor_aux(UIntType &result, UIntType &value) {
    constexpr std::size_t Size = std::min(Shift, End - Position);
    constexpr auto m = mask<UIntType>(Size, Position);
    UIntType const part = value & m;
    result |= part;

    constexpr auto NewShift = Position + Size;
    if constexpr (NewShift < End) {
        value ^= (part << Shift);
        unshift_left_xor_aux<End, Shift, NewShift>(result, value);
    }
}

} // namespace detail

template <std::size_t End, std::size_t Shift, class UIntType>
inline auto unshift_left_xor(UIntType value) -> UIntType {
    UIntType result = 0;
    detail::unshift_left_xor_aux<End, Shift>(result, value);
    return result;
}

template <std::size_t End, std::size_t Shift, std::size_t Position = 0,
          class UIntType>
inline auto unshift_left_xor(UIntType &result, UIntType value) {
    constexpr auto NewShift = Position + Shift;
    if constexpr (NewShift < End) {
        constexpr std::size_t Size = std::min(Shift, End - NewShift);
        constexpr auto m = mask<UIntType>(Size, Position);
        UIntType part = result & m;

        value ^= (part << Shift);
        detail::unshift_left_xor_aux<End, Shift, NewShift>(result, value);
    }
}

namespace detail {

template <std::size_t End, std::size_t Shift, std::size_t Position = 0,
          class UIntType>
inline auto unshift_right_xor_aux(UIntType &result, UIntType &value) {
    constexpr std::size_t Size = std::min(Shift, End - Position);
    constexpr auto m = mask<UIntType>(Size, End - Position - Size);
    UIntType const part = value & m;
    result |= part;

    constexpr auto NewShift = Position + Size;
    if constexpr (NewShift < End) {
        value ^= (part >> Shift);
        unshift_right_xor_aux<End, Shift, NewShift>(result, value);
    }
}

} // namespace detail

template <std::size_t End, std::size_t Shift, class UIntType>
inline auto unshift_right_xor(UIntType value) -> UIntType {
    UIntType result = 0;
    detail::unshift_right_xor_aux<End, Shift>(result, value);
    return result;
}

template <std::size_t End, int Shift, class UIntType>
inline auto unshift_xor(UIntType value) -> UIntType {
    if constexpr (Shift > 0) {
        return unshift_left_xor<End, Shift>(value);
    } else if constexpr (Shift < 0) {
        return unshift_right_xor<End, -Shift>(value);
    } else {
        return value; // No shift
    }
}

namespace detail {

template <std::size_t End, std::size_t Shift, std::size_t Position = 0,
          class UIntType>
inline auto unshift_left_plus_aux(UIntType &result, UIntType &value) {
    constexpr std::size_t Size = std::min(Shift, End - Position);
    constexpr auto m = mask<UIntType>(Size, Position);
    UIntType const part = value & m;
    result |= part;

    constexpr auto NewShift = Position + Size;
    if constexpr (NewShift < End) {
        value -= (part << Shift);
        unshift_left_plus_aux<End, Shift, NewShift>(result, value);
    }
}

} // namespace detail

template <std::size_t End, std::size_t Shift, class UIntType>
inline auto unshift_left_plus(UIntType value) -> UIntType {
    UIntType result = 0;
    detail::unshift_left_plus_aux<End, Shift>(result, value);
    return result;
}

} // namespace adhoc

#endif // BRNG_TOOLS_UNSHIFT
