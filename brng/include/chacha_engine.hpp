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

#include "tools/mask.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t R>
class chacha_engine final {
    static_assert(w <= std::numeric_limits<UIntType>::digits);
    static_assert(w <= 64);

    inline void generate_block() {
        constexpr std::array<std::uint32_t, 4> constants = {
            0x61707865, 0x3320646e, 0x79622d32, 0x6b206574};

        std::array<std::uint32_t, 16> input{};
        for (std::size_t i = 0; i < 4U; ++i) {
            input[i] = constants[i];
        }

        for (std::size_t i = 0; i < 8U; ++i) {
            input[4 + i] = this->keysetup[i];
        }

        input[12] = (this->ctr / 16U) & 0xffffffffU;
        input[13] = (this->ctr / 16U) >> 32U;
        input[14] = input[15] = 0xdeadbeefU;

        for (std::size_t i = 0; i < 16; ++i) {
            this->block[i] = input[i];
        }

        auto quarterroud =
            []<std::size_t a, std::size_t b, std::size_t c, std::size_t d>(
                std::array<std::uint32_t, 16> &x) {
                auto rotation = []<std::size_t ww, std::size_t n>(UIntType &x) {
                    x = (x << n) | (x >> (ww - n));
                };

                x[a] = x[a] + x[b];
                x[d] ^= x[a];
                rotation.template operator()<32, 16>(x[d]);
                x[c] = x[c] + x[d];
                x[b] ^= x[c];
                rotation.template operator()<32, 12>(x[b]);
                x[a] = x[a] + x[b];
                x[d] ^= x[a];
                rotation.template operator()<32, 8>(x[d]);
                x[c] = x[c] + x[d];
                x[b] ^= x[c];
                rotation.template operator()<32, 7>(x[b]);
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
        }
    }

  public:
    using result_type = UIntType;
    static constexpr std::uint64_t default_seed_1 = 0xb504f333f9de6484UL;
    static constexpr std::uint64_t default_seed_2 = 0UL;

    chacha_engine() : chacha_engine(default_seed_1, default_seed_2) {}

    chacha_engine(std::uint64_t seedval, std::uint64_t stream) {

        this->keysetup[0] = seedval & 0xffffffffU;
        this->keysetup[1] = seedval >> 32U;
        this->keysetup[2] = keysetup[3] = 0xdeadbeefU;
        this->keysetup[4] = stream & 0xffffffffU;
        this->keysetup[5] = stream >> 32U;
        this->keysetup[6] = keysetup[7] = 0xdeadbeefU;
        this->operator()<true>();
        this->operator()<false>();
    }

    template <class SeedSeq> explicit chacha_engine(SeedSeq &seq) {
        seq.generate(this->keysetup.begin(), this->keysetup.end());
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            std::uint64_t idx = this->ctr % 16U;
            if (idx == 0) {
                this->generate_block();
            }
            ++this->ctr;
            return this->block[idx];
        } else {
            --this->ctr;
            std::uint64_t idx = this->ctr % 16U;
            result_type result = this->block[idx];
            if (idx == 0) {
                this->ctr -= 16;
                this->generate_block();
                this->ctr += 16;
            }
            return result;
        }
    }

    template <bool FwdDirection = true> void discard(unsigned long long z) {
        std::uint64_t cycle = this->ctr / 16U;
        if constexpr (FwdDirection) {
            this->ctr += z;
        } else {
            this->ctr -= z;
        }
        std::uint64_t cycle_new = this->ctr / 16U;
        if (cycle != cycle_new) {
            --this->ctr;
            this->generate_block();
            ++this->ctr;
        }
    }

    static constexpr auto min() -> result_type {
        return static_cast<result_type>(0U);
    };

    static constexpr auto max() -> result_type { return mask<UIntType, w>(); };

    auto operator==(const chacha_engine &rhs) -> bool {
        return (this->keysetup == rhs.keysetup) && (this->ctr == rhs.ctr) &&
               (this->block == rhs.block);
    }

    auto operator!=(const chacha_engine &rhs) -> bool {
        return !(this->operator==(rhs));
    }

  private:
    std::array<std::uint32_t, 16> block{};
    std::array<std::uint32_t, 8> keysetup{};
    std::uint64_t ctr = 0;
};

using chacha4r = chacha_engine<std::uint32_t, 32, 4>;
using chacha5r = chacha_engine<std::uint32_t, 32, 5>;
using chacha6r = chacha_engine<std::uint32_t, 32, 6>;
using chacha8r = chacha_engine<std::uint32_t, 32, 8>;
using chacha20r = chacha_engine<std::uint32_t, 32, 20>;

} // namespace adhoc

#endif // BRNG_CHACHA_ENGINE
