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

#ifndef BRNG_TOOLS_MASK
#define BRNG_TOOLS_MASK

#include <cstddef>
#include <limits>
#include <type_traits>

namespace adhoc {

template <class UIntType>
constexpr auto mask(std::size_t Size, std::size_t Start = 0) -> UIntType {
// TODO: make this work on windows
#ifndef _MSC_VER
    static_assert(std::is_unsigned_v<UIntType>,
                  "UIntType must be an unsigned type");
#endif

    // Create a mask of `Size` bits set to 1, shifted to `Start` position
    return (static_cast<UIntType>(~UIntType{0}) >>
            (std::numeric_limits<UIntType>::digits - Size))
           << Start;
}
} // namespace adhoc

#endif // BRNG_TOOLS_MASK
