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

#ifndef COMBINATORICS_PARTITION_FUNCTION_HPP
#define COMBINATORICS_PARTITION_FUNCTION_HPP

namespace adhoc4 {

namespace detail {
constexpr auto partition_function_aux(int n, int k) -> int {
    if (k == 0) {
        return 0;
    }
    if (n == 0) {
        return 1;
    }
    if (n < 0) {
        return 0;
    }

    return partition_function_aux(n, k - 1) + partition_function_aux(n - k, k);
}
} // namespace detail

constexpr auto partition_function(int n) -> int {
    // convention
    if (n == 0) {
        return 1;
    }
    return detail::partition_function_aux(n, n);
}

} // namespace adhoc4

#endif // COMBINATORICS_PARTITION_FUNCTION_HPP
