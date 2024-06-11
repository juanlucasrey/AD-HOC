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

#ifndef ADHOC4_UTILS_POW_HPP
#define ADHOC4_UTILS_POW_HPP

#include <array>
#include <cmath>

namespace adhoc4::detail {

template <std::size_t N, class T> inline auto pow(T base) -> T {
    if constexpr (N == 1) {
        return base;
    } else if constexpr (N % 2 == 0) {
        T power = pow<N / 2>(base);
        return power * power;

    } else {
        T power = pow<N / 2>(base);
        return base * power * power;
    }
};

} // namespace adhoc4::detail

#endif // ADHOC4_UTILS_POW_HPP
