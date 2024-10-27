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

#ifndef ADHOC4_UTILS_ARRAY_CHAR_HPP
#define ADHOC4_UTILS_ARRAY_CHAR_HPP

#include <array>
#include <cstddef>

namespace adhoc4 {

namespace detail {

template <class IT, class Intpu1, class... Input>
constexpr void concatenate_aux(IT it, Intpu1 in, Input... ins) {
    it = std::copy_n(in.data(), in.size(), it);
    if constexpr (sizeof...(Input)) {
        concatenate_aux(it, ins...);
    }
}

} // namespace detail

template <class... Input> constexpr auto concatenate(Input... ins) {
    constexpr std::size_t N = (ins.size() + ...);
    std::array<char, N> result = {};
    detail::concatenate_aux(result.begin(), ins...);
    return result;
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_ARRAY_CHAR_HPP
