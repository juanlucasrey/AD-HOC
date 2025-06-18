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

#ifndef BRNG_TOOLS_BIT
#define BRNG_TOOLS_BIT

#include "mask.hpp"

#include <bit>
#include <limits>

namespace adhoc {

template <std::size_t w, class UIntType>
constexpr auto rotl(UIntType x, int s) noexcept -> UIntType {
    static_assert(w <= std::numeric_limits<UIntType>::digits);
    if constexpr (w == std::numeric_limits<UIntType>::digits) {
        return std::rotl(x, s);
    } else {
        constexpr UIntType m = mask<UIntType>(w);
        return ((x << s) & m) | (x >> (w - s));
    }
};

template <std::size_t w, class UIntType>
constexpr auto rotr(UIntType x, int s) noexcept -> UIntType {
    static_assert(w <= std::numeric_limits<UIntType>::digits);
    if constexpr (w == std::numeric_limits<UIntType>::digits) {
        return std::rotr(x, s);
    } else {
        constexpr UIntType m = mask<UIntType>(w);
        return (x >> s) | ((x << (w - s)) & m);
    }
};

template <class UIntType>
constexpr auto shift(UIntType x, int s) noexcept -> UIntType {
    if (s >= 0) {
        return x << s;
    } else {
        return x >> -s;
    }
};

} // namespace adhoc

#endif // BRNG_TOOLS_BIT
