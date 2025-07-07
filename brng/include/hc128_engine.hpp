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

#ifndef BRNG_HC128_ENGINE
#define BRNG_HC128_ENGINE

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
class hc128_engine final
    : public common_engine<UIntType, 32, hc128_engine<UIntType>> {
  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = 32U;
    static constexpr UIntType default_seed = 5489U;

  private:
    void init(std::array<UIntType, 8> const &init_array) {
        std::array<std::uint32_t, 1024 + 256> W{};

        for (std::size_t i = 0; i < 4; i++) {
            W[i] = ((std::uint32_t *)init_array.data())[i];
            W[i + 4] = ((std::uint32_t *)init_array.data())[i];
        }

        for (std::size_t i = 0; i < 4; i++) {
            W[i + 8] = ((std::uint32_t *)init_array.data())[i + 4];
            W[i + 12] = ((std::uint32_t *)init_array.data())[i + 4];
        }

        for (std::size_t i = 16; i < 1024 + 256; i++) {
            auto const tmp1 = std::rotr(W[i - 15], 7) ^
                              std::rotr(W[i - 15], 18) ^ (W[i - 15] >> 3U);

            auto const tmp2 = std::rotr(W[i - 2], 17) ^
                              std::rotr(W[i - 2], 19) ^ (W[i - 2] >> 10U);

            W[i] = tmp2 + W[i - 7] + tmp1 + W[i - 16] + i;
        }

        for (std::size_t i = 0; i < 512; i++, ++P) {
            P.at() = W[i + 256];
        }

        for (std::size_t i = 0; i < 512; i++, ++Q) {
            Q.at() = W[i + 256 + 512];
        }

        constexpr auto mask_idx = mask<std::size_t>(10);
        this->idx = mask_idx;
        --P;
        --Q;
        for (std::size_t i = 0; i < 1024; i++) {
            auto tmpidx = this->idx;
            ++tmpidx;
            tmpidx &= mask_idx;

            auto &val = tmpidx < 512U ? P.template at<1>() : Q.template at<1>();
            val = *(++(*this));
        }
        this->operator++();
    }

  public:
    hc128_engine() : hc128_engine(default_seed) {}

    explicit hc128_engine(UIntType s) {
        constexpr auto mask_res = hc128_engine::max();

        std::array<UIntType, 8> init_array{};
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

    template <class SeedSeq> explicit hc128_engine(SeedSeq &seq) {
        std::array<UIntType, 8> init_array{};
        seq.generate(init_array.begin(), init_array.end());
        this->init(init_array);
    }

    inline auto operator*() const -> value_type {
        constexpr auto mask_res = hc128_engine::max();
        bool const use_first = this->idx < 512U;
        auto const &PQ = use_first ? P : Q;
        auto const &QP_data = use_first ? this->Q.data() : this->P.data();
        auto const u = PQ.template at<512 - 12>();
        auto const a = static_cast<std::uint8_t>(u);
        auto const c = static_cast<std::uint8_t>(u >> 16U);
        UIntType result = (QP_data[a] + QP_data[256 + c]) ^ PQ.at();
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            result &= mask_res;
        }
        return result;
    }

    inline auto operator++() -> hc128_engine & {
        constexpr auto mask_idx = mask<std::size_t>(10);
        constexpr auto mask_res = hc128_engine::max();

        ++this->idx;
        this->idx &= mask_idx;
        bool const use_first = this->idx < 512U;
        auto &PQ = use_first ? P : Q;
        ++PQ;

        int const rot1 = use_first ? 10 : -10;
        int const rot2 = use_first ? 23 : -23;
        int const rot3 = use_first ? 8 : -8;
        PQ.at() = PQ.at() +
                  (rot<word_size>(PQ.template at<512 - 3>(), rot1) ^
                   rot<word_size>(PQ.template at<1>(), rot2)) +
                  rot<word_size>(PQ.template at<512 - 10>(), rot3);

        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            PQ.at() &= mask_res;
        }

        return *this;
    }

    inline auto operator--() -> hc128_engine & {
        constexpr auto mask_idx = mask<std::size_t>(10);
        constexpr auto mask_res = hc128_engine::max();
        bool const use_first = this->idx < 512U;
        auto &PQ = use_first ? P : Q;
        int const rot1 = use_first ? 10 : -10;
        int const rot2 = use_first ? 23 : -23;
        int const rot3 = use_first ? 8 : -8;
        PQ.at() = PQ.at() -
                  (rot<word_size>((PQ.template at<512 - 3>()), rot1) ^
                   rot<word_size>((PQ.template at<1>()), rot2)) -
                  rot<word_size>((PQ.template at<512 - 10>()), rot3);
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            PQ.at() &= mask_res;
        }
        --PQ;
        --this->idx;
        this->idx &= mask_idx;
        return *this;
    }

    using common_engine<UIntType, 32, hc128_engine>::operator++;
    using common_engine<UIntType, 32, hc128_engine>::operator--;
    using common_engine<UIntType, 32, hc128_engine>::operator==;

    auto operator==(const hc128_engine &rhs) const -> bool {
        return (this->idx == rhs.idx) && (this->P == rhs.P) &&
               (this->Q == rhs.Q);
    }

  private:
    std::size_t idx{0};
    circular_buffer<UIntType, 512> P{};
    circular_buffer<UIntType, 512> Q{};
};

} // namespace adhoc

#endif // BRNG_HC256_ENGINE
