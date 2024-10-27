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

#ifndef ADHOC4_UTILS_UNIVARIATE_HPP
#define ADHOC4_UTILS_UNIVARIATE_HPP

#include "../combinatorics/factorial.hpp"
#include "../differential_operator.hpp"
#include "../utils/index_sequence.hpp"

#include <array>
#include <cstddef>

namespace adhoc4 {

namespace detail {
template <std::size_t ToPow, std::size_t N, std::size_t N2 = N,
          std::size_t MaxOrder>
inline auto elevate_univariate_pow_2_aux(std::array<double, MaxOrder> const &ua,
                                         std::array<double, MaxOrder> &ua_res) {
    constexpr auto Pos = N2 + N + 1;
    constexpr auto EarlyStop = ToPow - 2;
    if constexpr (Pos < (MaxOrder - EarlyStop)) {
        if constexpr (N == N2) {
            ua_res[Pos] += ua[N] * ua[N2];
        } else {
            ua_res[Pos] += 2 * ua[N] * ua[N2];
        }
        elevate_univariate_pow_2_aux<ToPow, N, N2 + 1>(ua, ua_res);
    }
}

template <std::size_t ToPow, std::size_t N = 0, std::size_t MaxOrder>
inline auto elevate_univariate_pow_2(std::array<double, MaxOrder> const &ua,
                                     std::array<double, MaxOrder> &ua_res) {
    if constexpr (N < (MaxOrder / 2)) {
        elevate_univariate_pow_2_aux<ToPow, N>(ua, ua_res);
        elevate_univariate_pow_2<ToPow, N + 1>(ua, ua_res);
    }
}

template <std::size_t PrevPow, std::size_t ToPow, std::size_t N,
          std::size_t Loc = PrevPow, std::size_t MaxOrder>
inline auto
elevate_univariate_pow_3_up_aux(std::array<double, MaxOrder> const &ua,
                                std::array<double, MaxOrder> const &ua_prev,
                                std::array<double, MaxOrder> &ua_res) {
    constexpr auto EarlyStop = ToPow - PrevPow - 1;
    constexpr auto Pos = Loc + N;
    if constexpr (Pos < (MaxOrder - EarlyStop)) {
        ua_res[Pos] += ua[N] * ua_prev[Loc - 1];
        elevate_univariate_pow_3_up_aux<PrevPow, ToPow, N, Loc + 1>(ua, ua_prev,
                                                                    ua_res);
    }
}

template <std::size_t PrevPow, std::size_t ToPow, std::size_t N = 0,
          std::size_t MaxOrder>
inline auto
elevate_univariate_pow_3_up(std::array<double, MaxOrder> const &ua,
                            std::array<double, MaxOrder> const &ua_prev,
                            std::array<double, MaxOrder> &ua_res) {
    if constexpr (N + PrevPow - 1 < MaxOrder) {
        elevate_univariate_pow_3_up_aux<PrevPow, ToPow, N>(ua, ua_prev, ua_res);
        elevate_univariate_pow_3_up<PrevPow, ToPow, N + 1>(ua, ua_prev, ua_res);
    }
}
} // namespace detail

template <std::size_t FromPow, std::size_t ToPow, std::size_t MaxOrder>
inline void elevate_univariate(std::array<double, MaxOrder> const &ua,
                               std::array<double, MaxOrder> &ua_res) {
    if constexpr (FromPow < ToPow) {
        if constexpr (FromPow == 0) {
            ua_res = ua;
        } else if constexpr (FromPow == 1) {
            ua_res.fill(0);
            detail::elevate_univariate_pow_2<ToPow>(ua, ua_res);
        } else {
            std::array<double, MaxOrder> temp{0};
            detail::elevate_univariate_pow_3_up<FromPow, ToPow>(ua, ua_res,
                                                                temp);
            ua_res = temp;
        }

        elevate_univariate<FromPow + 1, ToPow>(ua, ua_res);
    }
}

namespace detail {
template <std::size_t Idx = 1, std::size_t MaxOrder>
inline void derivative_to_taylor_aux(std::array<double, MaxOrder> &ua) {
    if constexpr (Idx < MaxOrder) {
        constexpr double inverse_factorial = 1. / factorial(Idx + 1);
        ua[Idx] *= inverse_factorial;
        derivative_to_taylor_aux<Idx + 1>(ua);
    }
}
} // namespace detail

template <std::size_t MaxOrder>
inline void derivative_to_taylor(std::array<double, MaxOrder> &ua) {
    if constexpr (1 < MaxOrder) {
        detail::derivative_to_taylor_aux(ua);
    }
}

namespace detail {
template <class Id, std::size_t... I>
constexpr auto expand_univariate_aux(std::index_sequence<I...> /* i */) {
    // return std::tuple<std::tuple<der::d<I, Id>>...>{};
    return std::tuple<const der::d<I, Id>...>{};
}
} // namespace detail

template <class Id, std::size_t From, std::size_t To>
constexpr auto expand_univariate() {
    return detail::expand_univariate_aux<Id>(make_index_sequence<From, To>());
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_UNIVARIATE_HPP
