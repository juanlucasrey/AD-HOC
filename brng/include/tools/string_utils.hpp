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

#ifndef BRNG_TOOLS_STRING_UTILS
#define BRNG_TOOLS_STRING_UTILS

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace adhoc {

template <std::size_t Size>
constexpr auto constants_from_phrase(std::string_view str) {
    std::array<char, Size> arr;
    std::copy_n(str.data(), str.size(), arr.data());

    constexpr std::size_t NewSize = Size / 4;
    std::array<std::uint32_t, NewSize> res;
    for (std::size_t i = 0; i < NewSize; ++i) {
        res[i] = (static_cast<std::uint32_t>(arr[i * 4])) +
                 (static_cast<std::uint32_t>(arr[(i * 4) + 1]) << 8U) +
                 (static_cast<std::uint32_t>(arr[(i * 4) + 2]) << 16U) +
                 (static_cast<std::uint32_t>(arr[(i * 4) + 3]) << 24U);
    }
    return res;
}

} // namespace adhoc

#endif // BRNG_TOOLS_STRING_UTILS
