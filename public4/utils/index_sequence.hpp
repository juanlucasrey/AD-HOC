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

#ifndef ADHOC4_UTILS_INDEX_SEQUENCE_HPP
#define ADHOC4_UTILS_INDEX_SEQUENCE_HPP

#include <array>
#include <cmath>

namespace adhoc4 {

template <std::size_t Idx, class T, T... I>
constexpr auto get(std::integer_sequence<T, I...> /* i */) {
    return std::array<T, sizeof...(I)>{I...}[Idx];
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_INDEX_SEQUENCE_HPP
