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

#ifndef BRNG_HC256_ENGINE
#define BRNG_HC256_ENGINE

#include "tools/bit.hpp"
#include "tools/circular_buffer.hpp"
#include "tools/common_engine.hpp"
#include "tools/mask.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType>
class hc256_engine final
    : public common_engine<UIntType, 32, hc256_engine<UIntType>> {
  public:
    using value_type = UIntType;

    static constexpr std::size_t word_size = 32U;
    static constexpr UIntType default_seed = 5489U;

  private:
    void init(std::array<UIntType, 16> const &init_array) {
        std::array<std::uint32_t, 2048 + 512> W{};

        for (std::size_t i = 0; i < 16; i++) {
            W[i] = ((std::uint32_t *)init_array.data())[i];
        }

        for (std::size_t i = 16; i < 2048 + 512; i++) {
            auto const tmp1 = std::rotr(W[i - 15], 7) ^
                              std::rotr(W[i - 15], 18) ^ (W[i - 15] >> 3U);

            auto const tmp2 = std::rotr(W[i - 2], 17) ^
                              std::rotr(W[i - 2], 19) ^ (W[i - 2] >> 10U);

            W[i] = tmp2 + W[i - 7] + tmp1 + W[i - 16] + i;
        }

        for (std::size_t i = 0; i < 1024; i++, ++P) {
            P.at() = W[i + 512];
        }

        for (std::size_t i = 0; i < 1024; i++, ++Q) {
            Q.at() = W[i + 512 + 1024];
        }

        constexpr auto mask_idx = mask<std::size_t>(11);
        this->idx = mask_idx;
        --P;
        --Q;
        for (std::size_t i = 0; i < 4096; i++) {
            this->operator++();
        }
        this->operator++();
    }

  public:
    hc256_engine() : hc256_engine(default_seed) {}

    explicit hc256_engine(UIntType s) {
        constexpr auto mask_res = hc256_engine::max();

        std::array<UIntType, 16> init_array{};
        init_array[0] = UIntType(s);
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            init_array[0] &= mask_res;
        }

        if constexpr (word_size < std::numeric_limits<UIntType>::digits) {
            init_array[1] = UIntType(s >> 32);
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                init_array[1] &= mask_res;
            }
        }

        this->init(init_array);
    }

    template <class SeedSeq> explicit hc256_engine(SeedSeq &seq) {
        std::array<UIntType, 16> init_array{};
        seq.generate(init_array.begin(), init_array.end());
        this->init(init_array);
    }

    inline auto operator*() const -> value_type {
        constexpr auto mask_res = hc256_engine::max();
        bool const use_first = this->idx < 1024U;
        auto const &PQ = use_first ? P : Q;
        auto const &QP_data = use_first ? this->Q.data() : this->P.data();
        auto const u = PQ.template at<1024 - 12>();
        auto const a = static_cast<std::uint8_t>(u);
        auto const b = static_cast<std::uint8_t>(u >> 8U);
        auto const c = static_cast<std::uint8_t>(u >> 16U);
        auto const d = static_cast<std::uint8_t>(u >> 24U);
        std::uint32_t result = (QP_data[a] + QP_data[256 + b] +
                                QP_data[512 + c] + QP_data[768 + d]) ^
                               PQ.at();
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            result &= mask_res;
        }
        return result;
    }

    inline auto operator++() -> hc256_engine & {
        constexpr auto mask_idx = mask<std::size_t>(11);
        constexpr auto mask_half_idx = mask<std::size_t>(10);
        constexpr auto mask_res = hc256_engine::max();

        ++this->idx;
        this->idx &= mask_idx;
        bool const use_first = this->idx < 1024U;
        auto &PQ = use_first ? P : Q;
        ++PQ;

        auto const &QP_data = use_first ? this->Q.data() : this->P.data();
        PQ.at() +=
            PQ.template at<1024 - 10>() +
            (rot<word_size>(PQ.template at<1024 - 3>(), 10) ^
             rot<word_size>(PQ.template at<1>(), 23)) +
            QP_data[(PQ.template at<1024 - 3>() ^ (PQ.template at<1>())) &
                    mask_half_idx];

        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            PQ.at() &= mask_res;
        }

        return *this;
    }

    inline auto operator--() -> hc256_engine & {
        constexpr auto mask_idx = mask<std::size_t>(11);
        constexpr auto mask_half_idx = mask<std::size_t>(10);
        constexpr auto mask_res = hc256_engine::max();
        bool const use_first = this->idx < 1024U;
        auto &PQ = use_first ? P : Q;
        auto const &QP_data = use_first ? this->Q.data() : this->P.data();
        PQ.at() -=
            PQ.template at<1024 - 10>() +
            (rot<word_size>(PQ.template at<1024 - 3>(), 10) ^
             rot<word_size>(PQ.template at<1>(), 23)) +
            QP_data[(PQ.template at<1024 - 3>() ^ (PQ.template at<1>())) &
                    mask_half_idx];
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            PQ.at() &= mask_res;
        }
        --PQ;
        --this->idx;
        this->idx &= mask_idx;
        return *this;
    }

    using common_engine<UIntType, 32, hc256_engine>::operator++;
    using common_engine<UIntType, 32, hc256_engine>::operator--;
    using common_engine<UIntType, 32, hc256_engine>::operator==;

    auto operator==(const hc256_engine &rhs) const -> bool {
        return (this->idx == rhs.idx) && (this->P == rhs.P) &&
               (this->Q == rhs.Q);
    }

  private:
    std::size_t idx{0};
    circular_buffer<UIntType, 1024> P{};
    circular_buffer<UIntType, 1024> Q{};
};

} // namespace adhoc

#endif // BRNG_HC256_ENGINE
