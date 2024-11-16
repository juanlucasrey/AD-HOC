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

#ifndef COMBINATORICS_MULTINOMIAL_COEFFICIENT_HPP
#define COMBINATORICS_MULTINOMIAL_COEFFICIENT_HPP

#include <array>
#include <cstddef>

namespace adhoc4 {

template <std::size_t Bins, std::size_t Balls>
constexpr auto FirstMultinomial() {
    std::array<std::size_t, Bins> arr{0};
    arr.front() = Balls;
    return arr;
}

template <std::size_t Bins, std::size_t Balls>
constexpr auto LastMultinomial() {
    std::array<std::size_t, Bins> arr{0};
    arr.back() = Balls;
    return arr;
}

template <std::size_t Bins>
constexpr auto NextMultinomial(std::array<std::size_t, Bins> const &prev) {
    // check the previous mutlinomial coefficient is not the last
    std::array<std::size_t, Bins> arrlast{0};
    arrlast.back() = prev.back();
    if (arrlast == prev) {
        return prev;
    }

    std::array<std::size_t, Bins> arr = prev;

    std::size_t rem_val = arr.back();
    arr.back() = 0;
    std::size_t k = Bins - 2;
    while (arr[k] == 0 && k > 0) {
        k--;
    }
    arr[k]--;
    rem_val++;
    arr[k + 1] = rem_val;
    return arr;
}

template <std::size_t Order>
constexpr auto factorialarr() -> std::array<std::size_t, Order + 1> {
    std::array<std::size_t, Order + 1> res{};
    res[0] = 1;

    std::size_t k = 1;
    while (k <= Order) {
        res[k] = res[k - 1] * k;
        k++;
    }

    return res;
}

template <std::size_t Order, std::size_t Bins, std::size_t FactorialInputs>
constexpr auto
MultinomialCoeff(std::array<std::size_t, FactorialInputs> const &factorials,
                 std::array<std::size_t, Bins> const &coeffs) -> std::size_t {
    std::size_t res = factorials[Order];

    std::size_t k = 0;
    while (k < Bins) {
        res /= factorials[coeffs[k]];
        k++;
    }
    return res;
}

} // namespace adhoc4

#endif // COMBINATORICS_MULTINOMIAL_COEFFICIENT_HPP
