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

#ifndef BRNG_LFSR_XOR_ENGINE
#define BRNG_LFSR_XOR_ENGINE

#include "tools/linear_feedback_shift_register.hpp"
#include "tools/mask.hpp"

#include <array>
#include <cstdint>
#include <tuple>

namespace adhoc {

namespace detail {

template <std::size_t... Ints> struct make_tuple_of_triplets;

template <std::size_t F> struct make_tuple_of_triplets<F> {
    using type = std::tuple<>;
};

template <std::size_t F, std::size_t A, std::size_t B, std::size_t C,
          std::size_t... Rest>
struct make_tuple_of_triplets<F, A, B, C, Rest...> {
    using type = decltype(std::tuple_cat(
        std::tuple<std::index_sequence<F, A, B, C>>{},
        typename make_tuple_of_triplets<F, Rest...>::type{}));
};

template <std::size_t... Ints>
using make_tuple_of_triplets_t = typename make_tuple_of_triplets<Ints...>::type;

template <class... Ints> struct make_tuple_of_lfsr;

template <class UIntType> struct make_tuple_of_lfsr<UIntType, std::tuple<>> {
    using type = std::tuple<>;
};

template <class UIntType, std::size_t F, std::size_t A, std::size_t B,
          std::size_t C, class... Rest>
struct make_tuple_of_lfsr<
    UIntType, std::tuple<std::index_sequence<F, A, B, C>, Rest...>> {
    using type = decltype(std::tuple_cat(
        std::tuple<
            adhoc::linear_feedback_shift_register<UIntType, F, A, B, C>>{},
        typename make_tuple_of_lfsr<UIntType, std::tuple<Rest...>>::type{}));
};

template <class UIntType, std::size_t... Ints>
using make_tuple_of_lfsr_t =
    typename make_tuple_of_lfsr<UIntType,
                                make_tuple_of_triplets_t<Ints...>>::type;

} // namespace detail

template <class UIntType, std::size_t w, std::size_t... vals> class xor_engine {
    static_assert(sizeof...(vals) % 3 == 0);
    static_assert(w <= std::numeric_limits<UIntType>::digits);
    static_assert(std::is_unsigned_v<UIntType>);
    static constexpr UIntType number_of_lfsr = sizeof...(vals) / 3;

    template <std::size_t Idx = 0>
    void init_consistent(std::array<UIntType, number_of_lfsr> const &init) {
        std::get<Idx>(this->state) =
            std::tuple_element_t<Idx, tuple_lfsr>(std::get<Idx>(init));
        if constexpr (static_cast<UIntType>(Idx) < (number_of_lfsr - 1)) {
            this->init_consistent<Idx + 1>(init);
        }
    }

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr UIntType default_seed =
        w >= 64 ? 123456789123456789ULL : 987654321U;

    xor_engine() : xor_engine(default_seed) {}

    explicit xor_engine(result_type value) {

        if (value == 0) {
            value = default_seed;
        }

        std::array<UIntType, number_of_lfsr> init{};
        init.fill(value);
        this->init_consistent(init);
    };

    template <class SeedSeq> explicit xor_engine(SeedSeq &seq) {
        std::array<std::uint32_t, number_of_lfsr> generated_sequence;
        seq.generate(generated_sequence.begin(), generated_sequence.end());

        std::array<UIntType, number_of_lfsr> init;
        std::transform(
            generated_sequence.begin(), generated_sequence.end(), init.begin(),
            [](std::uint32_t v) { return static_cast<UIntType>(v); });

        this->init_consistent(init);
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            return std::apply(
                []<class... Args>(Args &...args) { return (args() ^ ...); },
                this->state);
        } else {
            return std::apply(
                []<class... Args>(Args &...args) {
                    return (args.template operator()<false>() ^ ...);
                },
                this->state);
        }
    }

    static constexpr auto min() -> UIntType {
        return static_cast<UIntType>(0U);
    }

    static constexpr auto max() -> UIntType {
        return mask<UIntType>(word_size);
    }

    auto operator==(const xor_engine &rhs) const -> bool {
        return this->state == rhs.state;
    }

    auto operator!=(const xor_engine &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    using tuple_lfsr = detail::make_tuple_of_lfsr_t<UIntType, w, vals...>;
    tuple_lfsr state;
};

template <class UIntType>
using taus88 = xor_engine<UIntType, 32, 31, 13, 12, 29, 2, 4, 28, 3, 17>;

template <class UIntType>
using lfsr113 =
    xor_engine<UIntType, 32, 31, 6, 18, 29, 2, 2, 28, 13, 7, 25, 3, 13>;

template <class UIntType>
using lfsr258 = xor_engine<UIntType, 64, 63, 1, 10, 55, 24, 5, 52, 3, 29, 47, 5,
                           23, 41, 3, 8>;

} // namespace adhoc

#endif // BRNG_LFSR_XOR_ENGINE
