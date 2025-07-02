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

#ifndef BRNG_DSFMT_ENGINE
#define BRNG_DSFMT_ENGINE

#include "tools/bit.hpp"
#include "tools/circular_buffer.hpp"
#include "tools/common_engine.hpp"
#include "tools/mask.hpp"
#include "tools/seed_seq_filler.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace adhoc {

template <class UIntType, std::size_t mexp, std::size_t pos1, std::size_t sl,
          std::uint64_t msk1, std::uint64_t msk2, std::uint64_t fix1,
          std::uint64_t fix2, std::uint64_t parity1, std::uint64_t parity2>
class dsfmt_engine final
    : public common_engine<UIntType, 52,
                           dsfmt_engine<UIntType, mexp, pos1, sl, msk1, msk2,
                                        fix1, fix2, parity1, parity2>> {
  private:
    static constexpr std::size_t n = ((mexp - 128) / 104) + 1;

    using w128_t = std::array<UIntType, 2>;

    void initial_mask() {
        constexpr UIntType low_mask = mask<UIntType>(52);
        constexpr UIntType high_mask = mask<UIntType>(10, 52);

        for (auto &word : this->state.data()) {
            for (auto &subword : word) {
                subword = (subword & low_mask) | high_mask;
            }
        }
    }

    void period_certification() {
        constexpr std::array<std::uint64_t, 2> parity = {parity1, parity2};
        std::uint64_t inner = ((this->lung[0] ^ fix1) & parity[0]) ^
                              ((this->lung[1] ^ fix2) & parity[1]);

        for (std::size_t i = 32; i > 0; i >>= 1U) {
            inner ^= inner >> i;
        }
        inner &= 1U;
        // check OK
        if (inner == 1U) {
            return;
        }

        // check NG, and modification
        if constexpr ((parity2 & 1U) == 1) {
            this->lung[1] ^= 1;
        } else {
            for (std::size_t i = 2; i-- > 0;) {
                std::uint64_t work = 1;
                for (std::size_t j = 0; j < 64; j++) {
                    if ((work & parity[i]) != 0) {
                        this->lung[i] ^= work;
                        return;
                    }
                    work <<= 1U;
                }
            }
        }
    }

    void init_consistent() {
        this->initial_mask();
        this->period_certification();
        --this->state;
        this->operator++();
    }

  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = 52;
    static constexpr UIntType default_seed = 5489U;

    dsfmt_engine() : dsfmt_engine(default_seed) {}

    explicit dsfmt_engine(std::uint32_t seed) {

        std::size_t i = 0;
        for (auto &arr : this->state.data()) {
            for (auto &val : arr) {
                val = static_cast<UIntType>(seed);
                seed = 1812433253UL * (seed ^ (seed >> 30U)) + ++i;

                val |= static_cast<UIntType>(seed) << 32U;
                seed = 1812433253UL * (seed ^ (seed >> 30U)) + ++i;
            }
        }

        for (auto &val : lung) {
            val = static_cast<UIntType>(seed);
            seed = 1812433253UL * (seed ^ (seed >> 30U)) + ++i;

            val |= static_cast<UIntType>(seed) << 32U;
            seed = 1812433253UL * (seed ^ (seed >> 30U)) + ++i;
        }

        this->init_consistent();
    }

    template <std::size_t key_length>
    explicit dsfmt_engine(std::array<std::uint32_t, key_length> init_key) {
        auto &state_data = this->state.data();
        constexpr std::size_t size = (n + 1) * 4;
        constexpr std::size_t lag = size >= 623  ? 11
                                    : size >= 68 ? 7
                                    : size >= 39 ? 5
                                                 : 3;
        constexpr std::size_t mid = (size - lag) / 2;

        for (auto &arr : state_data) {
            for (auto &val : arr) {
                val = 10055284024492657547UL;
            }
        }

        for (auto &val : lung) {
            val = 10055284024492657547UL;
        }

        constexpr std::size_t count = std::max(size, key_length + 1U);

        auto get_val32 = [&state_data, this](std::size_t i) -> std::uint32_t {
            std::size_t const state_idx = i / 4;
            std::size_t const sub_idx = i % 4;

            std::size_t const sub_idx_64 = sub_idx / 2;
            std::size_t const sub_idx_64_half = sub_idx % 2;
            std::uint64_t val = state_idx == state_data.size()
                                    ? this->lung[sub_idx_64]
                                    : state_data.data()[state_idx][sub_idx_64];
            if (sub_idx_64_half == 1) {
                val >>= 32U;
            }
            return static_cast<std::uint32_t>(val);
        };

        auto set_val32 = [&state_data, this](std::size_t i, std::uint32_t r) {
            std::size_t const state_idx = i / 4;
            std::size_t const sub_idx = i % 4;

            std::size_t const sub_idx_64 = sub_idx / 2;
            std::size_t const sub_idx_64_half = sub_idx % 2;
            std::uint64_t &val = state_idx == state_data.size()
                                     ? lung[sub_idx_64]
                                     : state_data.data()[state_idx][sub_idx_64];
            auto r64 = static_cast<std::uint64_t>(r);
            if (sub_idx_64_half == 1) {
                constexpr auto mask_lower_half = mask<std::uint64_t>(32);
                r64 <<= 32U;
                val = r64 | (val & mask_lower_half);
            } else {
                constexpr auto mask_upper_half = mask<std::uint64_t>(32, 32);
                val = r64 | (val & mask_upper_half);
            }
        };

        for (std::size_t j = 0; j < count; ++j) {
            std::size_t const i = j % size;
            std::uint32_t r = get_val32(i) ^ get_val32((i + mid) % size) ^
                              get_val32((i + size - 1) % size);

            r = (r ^ (r >> 27U)) * 1664525UL;
            set_val32((i + mid) % size, get_val32((i + mid) % size) + r);

            if (j == 0) {
                r += key_length;
            } else {
                r += i;
                if ((j - 1) < key_length) {
                    r += init_key[j - 1];
                }
            }

            set_val32((i + mid + lag) % size,
                      get_val32((i + mid + lag) % size) + r);
            set_val32(i, r);
        }

        for (std::size_t i = 0; i < size; i++) {
            std::uint32_t r = get_val32(i) + get_val32((i + mid) % size) +
                              get_val32((i + size - 1) % size);
            r = (r ^ (r >> 27U)) * 1566083941UL;

            set_val32((i + mid) % size, get_val32((i + mid) % size) ^ r);
            r -= i;
            set_val32((i + mid + lag) % size,
                      get_val32((i + mid + lag) % size) ^ r);
            set_val32(i, r);
        }
        this->init_consistent();
    }

    template <class SeedSeq> explicit dsfmt_engine(SeedSeq &seq) {

        seed_seq_filler<64, (n + 1) * 4> seq_filler(seq);

        for (auto &word : this->state.data()) {
            for (auto &subword : word) {
                subword = seq_filler.template get<UIntType>();
            }
        }

        for (auto &subword : this->lung) {
            subword = seq_filler.template get<UIntType>();
        }

        this->init_consistent();
    }

    inline auto operator*() const -> value_type {
        constexpr auto mask_result = dsfmt_engine::max();
        return this->state.at()[this->idx] & mask_result;
    };

    inline auto operator++() -> dsfmt_engine & {
        ++this->idx;
        if (this->idx == 2) {
            ++this->state;
            this->idx = 0;

            auto &current = this->state.at();
            auto const &word_pos1 = this->state.template at<pos1>();

            auto const L0_prv = this->lung[0];

            if constexpr (64 < std::numeric_limits<UIntType>::digits) {
                auto constexpr mask64 = mask<UIntType>(64);
                this->lung[0] = ((current[0] << sl) & mask64) ^
                                rotr<64>(this->lung[1], 32) ^ word_pos1[0];
                this->lung[1] = ((current[1] << sl) & mask64) ^
                                rotr<64>(L0_prv, 32) ^ word_pos1[1];
            } else {
                this->lung[0] = (current[0] << sl) ^
                                rotr<64>(this->lung[1], 32) ^ word_pos1[0];
                this->lung[1] =
                    (current[1] << sl) ^ rotr<64>(L0_prv, 32) ^ word_pos1[1];
            }

            static constexpr std::size_t DSFMT_SR = 64U - word_size;
            current[0] ^= (this->lung[0] >> DSFMT_SR) ^ (this->lung[0] & msk1);
            current[1] ^= (this->lung[1] >> DSFMT_SR) ^ (this->lung[1] & msk2);
        }
        return *this;
    }

    inline auto operator--() -> dsfmt_engine & {
        if (this->idx == 0) {
            auto &current = this->state.at();

            static constexpr std::size_t DSFMT_SR = 64U - word_size;
            current[0] ^= (this->lung[0] >> DSFMT_SR) ^ (this->lung[0] & msk1);
            current[1] ^= (this->lung[1] >> DSFMT_SR) ^ (this->lung[1] & msk2);

            auto const &word_pos1 = this->state.template at<pos1>();
            auto const L0_prv = this->lung[0];

            auto constexpr mask64 = mask<UIntType>(64);
            if constexpr (64 < std::numeric_limits<UIntType>::digits) {
                this->lung[0] = ((current[1] << sl) & mask64) ^ this->lung[1] ^
                                word_pos1[1];
            } else {
                this->lung[0] =
                    (current[1] << sl) ^ this->lung[1] ^ word_pos1[1];
            }
            this->lung[0] = rotr<64>(this->lung[0], 32);

            if constexpr (64 < std::numeric_limits<UIntType>::digits) {
                this->lung[1] =
                    ((current[0] << sl) & mask64) ^ L0_prv ^ word_pos1[0];
            } else {
                this->lung[1] = (current[0] << sl) ^ L0_prv ^ word_pos1[0];
            }
            this->lung[1] = rotr<64>(this->lung[1], 32);

            this->idx = 2;
            --this->state;
        }
        --this->idx;
        return *this;
    }

    using common_engine<UIntType, 52, dsfmt_engine>::operator++;
    using common_engine<UIntType, 52, dsfmt_engine>::operator--;
    using common_engine<UIntType, 52, dsfmt_engine>::operator==;

    auto operator==(const dsfmt_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->lung == rhs.lung) &&
               (this->idx == rhs.idx);
    }

  private:
    circular_buffer<w128_t, n> state{};
    w128_t lung{};

    std::size_t idx{2 - 1};
};

using dsfmt521 = dsfmt_engine<std::uint64_t, 521, 3, 25, 0x000fbfefff77efffUL,
                              0x000ffeebfbdfbfdfUL, 0xcfb393d661638469UL,
                              0xc166867883ae2adbUL, 0xccaa588000000000UL,
                              0x0000000000000001UL>;
using dsfmt1279 =
    dsfmt_engine<std::uint64_t, 1279, 9, 19, 0x000efff7ffddffee,
                 0x000fbffffff77fff, 0xb66627623d1a31be, 0x04b6c51147b6109b,
                 0x7049f2da382a6aeb, 0xde4ca84a40000001>;

using dsfmt2203 =
    dsfmt_engine<std::uint64_t, 2203, 7, 19, 0x000fdffff5edbfff,
                 0x000f77fffffffbfe, 0xb14e907a39338485, 0xf98f0735c637ef90,
                 0x8000000000000000, 0x0000000000000001>;

using dsfmt4253 =
    dsfmt_engine<std::uint64_t, 4253, 19, 19, 0x0007b7fffef5feff,
                 0x000ffdffeffefbfc, 0x80901b5fd7a11c65, 0x5a63ff0e7cb0ba74,
                 0x1ad277be12000000, 0x0000000000000001>;

using dsfmt11213 =
    dsfmt_engine<std::uint64_t, 11213, 37, 19, 0x000ffffffdf7fffd,
                 0x000dfffffff6bfff, 0xd0ef7b7c75b06793, 0x9c50ff4caae0a641,
                 0x8234c51207c80000, 0x0000000000000001>;

using dsfmt19937 =
    dsfmt_engine<std::uint64_t, 19937, 117, 19, 0x000ffafffffffb3f,
                 0x000ffdfffc90fffd, 0x90014964b32f4329, 0x3b8d12ac548a7c7a,
                 0x3d84e1ac0dc82880, 0x0000000000000001>;

using dsfmt44497 =
    dsfmt_engine<std::uint64_t, 44497, 304, 19, 0x000ff6dfffffffef,
                 0x0007ffdddeefff6f, 0x75d910f235f6e10e, 0x7b32158aedc8e969,
                 0x4c3356b2a0000000, 0x0000000000000001>;

using dsfmt86243 =
    dsfmt_engine<std::uint64_t, 86243, 231, 13, 0x000ffedff6ffffdf,
                 0x000ffff7fdffff7e, 0x1d553e776b975e68, 0x648faadf1416bf91,
                 0x5f2cd03e2758a373, 0xc0b7eb8410000001>;

using dsfmt132049 =
    dsfmt_engine<std::uint64_t, 132049, 371, 23, 0x000fb9f4eff4bf77,
                 0x000fffffbfefff37, 0x4ce24c0e4e234f3b, 0x62612409b5665c2d,
                 0x181232889145d000, 0x0000000000000001>;

using dsfmt216091 =
    dsfmt_engine<std::uint64_t, 216091, 1890, 23, 0x000bf7df7fefcfff,
                 0x000e7ffffef737ff, 0xd7f95a04764c27d7, 0x6a483861810bebc2,
                 0x3af0a8f3d5600000, 0x0000000000000001>;

} // namespace adhoc

#endif // BRNG_DSFMT_ENGINE
