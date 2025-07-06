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

#ifndef BRNG_CHACHA_ENGINE
#define BRNG_CHACHA_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/mask.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t R>
class chacha_engine final
    : public common_engine<UIntType, 32, chacha_engine<UIntType, R>> {
    static constexpr std::size_t w = 32;
    static_assert(32 <= std::numeric_limits<UIntType>::digits);

    inline void increase_counter() {
        ++this->ctr[0];
        if (this->ctr[0] == 0) {
            ++this->ctr[1];
        }
    }

    inline void increase_counter(int amount) {
        this->ctr[0] += amount;
        if (this->ctr[0] < amount) {
            ++this->ctr[1];
        }
    }

    inline void decrease_counter() {
        if (this->ctr[0] == 0) {
            --this->ctr[1];
        }
        --this->ctr[0];
    }

    inline void decrease_counter(int amount) {
        if (this->ctr[0] < amount) {
            --this->ctr[1];
        }
        this->ctr[0] -= amount;
    }

    inline void generate_block() {
        constexpr std::array<UIntType, 4> constants = {0x61707865, 0x3320646e,
                                                       0x79622d32, 0x6b206574};

        std::array<UIntType, 16> input{};
        for (std::size_t i = 0; i < 4U; ++i) {
            input[i] = constants[i];
        }

        for (std::size_t i = 0; i < 8U; ++i) {
            input[4 + i] = this->keysetup[i];
        }

        input[12] = static_cast<UIntType>((this->ctr[0] / 16U) & mask_result);

        input[13] = ((static_cast<UIntType>(this->ctr[1]) % 16) << 28U) |
                    ((this->ctr[0] / 16) >> 32);
        input[14] = static_cast<UIntType>((this->ctr[1] / 16) & mask_result);
        input[15] = static_cast<UIntType>((this->ctr[1] / 16) >> 32U);

        for (std::size_t i = 0; i < 16; ++i) {
            this->block[i] = input[i];
        }

        auto quarterroud =
            []<std::size_t a, std::size_t b, std::size_t c, std::size_t d>(
                std::array<UIntType, 16> &x) {
                x[a] = x[a] + x[b];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[a] &= mask_result;
                }
                x[d] ^= x[a];
                x[d] = rotl<w>(x[d], 16);
                x[c] = x[c] + x[d];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[c] &= mask_result;
                }
                x[b] ^= x[c];
                x[b] = rotl<w>(x[b], 12);
                x[a] = x[a] + x[b];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[a] &= mask_result;
                }
                x[d] ^= x[a];
                x[d] = rotl<w>(x[d], 8);
                x[c] = x[c] + x[d];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[c] &= mask_result;
                }
                x[b] ^= x[c];
                x[b] = rotl<w>(x[b], 7);
            };

        for (std::size_t i = 1; i < R; i += 2) {
            quarterroud.template operator()<0, 4, 8, 12>(this->block);
            quarterroud.template operator()<1, 5, 9, 13>(this->block);
            quarterroud.template operator()<2, 6, 10, 14>(this->block);
            quarterroud.template operator()<3, 7, 11, 15>(this->block);
            quarterroud.template operator()<0, 5, 10, 15>(this->block);
            quarterroud.template operator()<1, 6, 11, 12>(this->block);
            quarterroud.template operator()<2, 7, 8, 13>(this->block);
            quarterroud.template operator()<3, 4, 9, 14>(this->block);
        }

        if (R & 1U) {
            quarterroud.template operator()<0, 5, 10, 15>(this->block);
            quarterroud.template operator()<1, 6, 11, 12>(this->block);
            quarterroud.template operator()<2, 7, 8, 13>(this->block);
            quarterroud.template operator()<3, 4, 9, 14>(this->block);
        }

        for (std::size_t i = 0; i < 16; ++i) {
            this->block[i] += input[i];
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->block[i] &= mask_result;
            }
        }
    }

  public:
    using value_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::uint64_t default_seed_1 = 0xb504f333f9de6484UL;
    static constexpr UIntType mask_result = chacha_engine::max();

    chacha_engine() : chacha_engine(default_seed_1) {}

    chacha_engine(std::uint64_t seedval1, std::uint64_t seedval2 = 0,
                  std::uint64_t stream1 = 0, std::uint64_t stream2 = 0) {

        this->keysetup[0] = seedval1 & mask_result;
        this->keysetup[1] = seedval1 >> 32U;
        this->keysetup[2] = seedval2 & mask_result;
        this->keysetup[3] = seedval2 >> 32U;
        this->keysetup[4] = stream1 & mask_result;
        this->keysetup[5] = stream1 >> 32U;
        this->keysetup[6] = stream2 & mask_result;
        this->keysetup[7] = stream2 >> 32U;
        this->generate_block();
    }

    template <class SeedSeq> explicit chacha_engine(SeedSeq &seq) {
        seq.generate(this->keysetup.begin(), this->keysetup.end());
        this->generate_block();
    }

    inline auto operator*() const -> value_type {
        constexpr auto mask4 = mask<std::size_t>(4);
        return this->block[this->ctr[0] & mask4];
    };

    inline auto operator++() -> chacha_engine & {
        constexpr auto mask4 = mask<std::size_t>(4);

        bool const gen = (this->ctr[0] & mask4) == 15;
        this->increase_counter();
        if (gen) {
            this->generate_block();
        }
        return *this;
    }

    inline auto operator--() -> chacha_engine & {
        constexpr auto mask4 = mask<std::size_t>(4);

        this->decrease_counter();
        if ((this->ctr[0] & mask4) == 15) {
            this->generate_block();
        }
        return *this;
    }

    using common_engine<UIntType, w, chacha_engine>::operator++;
    using common_engine<UIntType, w, chacha_engine>::operator--;
    using common_engine<UIntType, w, chacha_engine>::operator==;

    auto operator==(const chacha_engine &rhs) const -> bool {
        return (this->keysetup == rhs.keysetup) && (this->ctr == rhs.ctr) &&
               (this->block == rhs.block);
    }

  private:
    std::array<UIntType, 16> block{};
    std::array<UIntType, 8> keysetup{};
    std::array<std::size_t, 2> ctr{};
};

using chacha20 = chacha_engine<std::uint_fast32_t, 20>;

} // namespace adhoc

#endif // BRNG_CHACHA_ENGINE
