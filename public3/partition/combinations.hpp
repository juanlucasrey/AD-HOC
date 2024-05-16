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

#ifndef ADHOC3_PARTITION_COMBINATIONS_HPP
#define ADHOC3_PARTITION_COMBINATIONS_HPP

#include "binomial_coefficient.hpp"

#include <cstddef>

namespace adhoc3 {

constexpr inline auto combinations(std::size_t bins, std::size_t balls) noexcept
    -> std::size_t {
    return binomial_coefficient(bins + balls - 1, bins - 1);
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_COMBINATIONS_HPP
