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

#ifndef COMBINATORICS_INTEGER_PARTITION_HPP
#define COMBINATORICS_INTEGER_PARTITION_HPP

#include <array>
#include <cstddef>

namespace adhoc {

template <std::size_t N> constexpr auto FirstPartition() {
    std::array<std::size_t, N> arr{0};
    arr.back() = 1;
    return arr;
}

template <std::size_t N> constexpr auto LastPartition() {
    std::array<std::size_t, N> arr{0};
    arr.front() = N;
    return arr;
}

template <std::size_t N>
constexpr auto NextPartition(std::array<std::size_t, N> const &prev) {
    if (prev == LastPartition<N>()) {
        return prev;
    }

    std::array<std::size_t, N> arr = prev;

    // Find the rightmost non-one value in arr. Also, update the
    // rem_val so that we know how much value can be accommodated
    std::size_t rem_val = arr.front();
    arr.front() = 0;
    std::size_t k = 1;
    while (arr[k] == 0) {
        k++;
    }

    arr[k]--;
    rem_val += (k + 1);

    // If rem_val is more, then the sorted order is violated. Divide
    // rem_val in different values of size p[k] and copy these values at
    // different positions after p[k]
    if (rem_val > k) {
        // std::div will be constexpr in C++23 so we use a manual version in
        // the meantime
        std::size_t const div = rem_val / k;
        arr[k - 1] += div;
        rem_val -= k * div;
    }

    if (rem_val) {
        arr[rem_val - 1]++;
    }

    return arr;
}

template <std::size_t N>
constexpr auto PrevPartition(std::array<std::size_t, N> const &prev) {
    if (prev == FirstPartition<N>()) {
        return prev;
    }

    std::array<std::size_t, N> arr = prev;

    // Find the rightmost non-one value in arr. Also, update the
    // rem_val so that we know how much value can be accommodated
    std::size_t rem_val = arr.front();
    arr.front() = 0;
    std::size_t k = 1;
    while (k >= rem_val) {
        rem_val += arr[k] * (k + 1);
        arr[k] = 0;
        k++;
    }

    arr[k]++;
    rem_val -= (k + 1);

    arr.front() = rem_val;
    return arr;
}

namespace detail {

constexpr auto factorial(std::size_t n) -> std::size_t {
    return n <= 1UL ? 1UL : (n * factorial(n - 1UL));
}

constexpr auto pow(std::size_t x, std::size_t n) -> std::size_t {
    if (n == 0UL) {
        return 1UL;
    }

    if ((n & 1UL) == 0UL) {
        return pow(x * x, n / 2UL);
    }

    return x * pow(x * x, (n - 1UL) / 2UL);
}

template <std::size_t N, std::size_t I>
constexpr auto BellDenom(std::array<std::size_t, N> const &str) -> std::size_t {
    if constexpr (I == 0) {
        return 1UL;
    } else {
        return pow(factorial(I), str[I - 1]) * factorial(str[I - 1]) *
               BellDenom<N, I - 1>(str);
    }
}

} // namespace detail

template <std::size_t N>
constexpr auto BellCoeff(std::array<std::size_t, N> const &str) -> std::size_t {
    return detail::factorial(N) / detail::BellDenom<N, N>(str);
}

} // namespace adhoc

#endif // COMBINATORICS_INTEGER_PARTITION_HPP
