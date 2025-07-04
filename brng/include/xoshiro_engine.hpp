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

#ifndef BRNG_XOSHIRO_ENGINE
#define BRNG_XOSHIRO_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/seed_seq_filler.hpp"
#include "tools/unshift.hpp"

#include <array>
#include <cstddef>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t n, std::size_t a,
          std::size_t b, std::size_t orot = 0, UIntType mult1 = 0,
          UIntType mult2 = 0, std::size_t idx1 = 2, std::size_t idx2 = 2>
class xoshiro_engine final
    : public common_engine<UIntType, w,
                           xoshiro_engine<UIntType, w, n, a, b, orot, mult1,
                                          mult2, idx1, idx2>> {
  private:
    void init4(UIntType s0, UIntType s1, UIntType s2, UIntType s3) {
        this->state[0] = s0 & mask_result;
        this->state[1] = s1 & mask_result;
        this->state[2] = s2 & mask_result;
        this->state[3] = s3 & mask_result;
        this->state[3] = (this->state[0] || this->state[1] || this->state[2] ||
                          this->state[3])
                             ? this->state[3]
                             : 1U;
    }

    void init8(UIntType s0, UIntType s1, UIntType s2, UIntType s3, UIntType s4,
               UIntType s5, UIntType s6, UIntType s7) {
        this->state[0] = s0 & mask_result;
        this->state[1] = s1 & mask_result;
        this->state[2] = s2 & mask_result;
        this->state[3] = s3 & mask_result;
        this->state[4] = s4 & mask_result;
        this->state[5] = s5 & mask_result;
        this->state[6] = s6 & mask_result;
        this->state[7] = s7 & mask_result;
        this->state[7] = (this->state[0] || this->state[1] || this->state[2] ||
                          this->state[3] || this->state[4] || this->state[5] ||
                          this->state[6] || this->state[7])
                             ? this->state[7]
                             : 1U;
    }

  public:
    static_assert(n == 4 || n == 8);
    static_assert(a < w);
    static_assert(b < w);

    using value_type = UIntType;
    static constexpr UIntType mask_result = xoshiro_engine::max();
    static constexpr UIntType default_seed = 1;

    xoshiro_engine() : xoshiro_engine(default_seed) {}

    explicit xoshiro_engine(UIntType seed) {

        if constexpr (n == 4) {
            this->init4(
                0x01d353e5f3993bb0 + seed, 0x7b9c0df6cb193b20 * (seed + 1),
                0xfdfcaa91110765b6 - seed, 0x2d24cbe0ef44dcd2 * (seed - 1));
        } else {
            this->init8(
                0x1ced436497db2a59 + seed, 0x75474f85d8a6892c * (seed + 1),
                0xa0fef4b8094c9c86 - seed, 0x748fa1a9bb555169 * (seed - 1),
                0xd7a59a6d64e66858 + seed, 0xf03b7efdb73db601 * (seed + 1),
                0xfab342a99dd71962 - seed, 0x8a6921456faa6b54 * (seed - 1));
        }

        this->discard(16);
    }

    template <class SeedSeq> explicit xoshiro_engine(SeedSeq &seq) {
        seed_seq_filler<w, n> seq_filler(seq);
        seq_filler.template generate<UIntType>(this->state.begin(),
                                               this->state.end());
    }

    inline auto operator*() const -> value_type {
        UIntType result = 0;

        if constexpr (mult1 && orot) {
            result = this->state[1] * mult1;
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                result &= mask_result;
            }

            result = rotl<w>(result, orot) * mult2;
        } else {
            result = this->state[0] + this->state[idx1];

            if constexpr (orot) {
                if constexpr (w != std::numeric_limits<UIntType>::digits) {
                    result &= mask_result;
                }
                result = rotl<w>(result, orot);
                result += this->state[idx2];
            }
        }

        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            result &= mask_result;
        }
        return result;
    }

    inline auto operator++() -> xoshiro_engine & {
        if constexpr (n == 8) {
            UIntType t = this->state[1] << a;
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                t &= mask_result;
            }
            this->state[2] ^= this->state[0];
            this->state[5] ^= this->state[1];
            this->state[1] ^= this->state[2];
            this->state[7] ^= this->state[3];
            this->state[3] ^= this->state[4];
            this->state[4] ^= this->state[5];
            this->state[0] ^= this->state[6];
            this->state[6] ^= this->state[7];
            this->state[6] ^= t;
            this->state[7] = rotl<w>(this->state[7], b);
        } else {
            UIntType t = this->state[1] << a;
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                t &= mask_result;
            }

            this->state[2] ^= this->state[0];
            this->state[3] ^= this->state[1];
            this->state[1] ^= this->state[2];
            this->state[0] ^= this->state[3];
            this->state[2] ^= t;
            this->state[3] = rotl<w>(this->state[3], b);
        }

        return *this;
    }

    inline auto operator--() -> xoshiro_engine & {
        if constexpr (n == 8) {
            this->state[7] = rotr<w>(this->state[7], b);
            this->state[1] ^= this->state[2];
            UIntType t = this->state[1] << a;
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                t &= mask_result;
            }
            this->state[6] ^= t;
            this->state[6] ^= this->state[7];
            this->state[0] ^= this->state[6];
            this->state[4] ^= this->state[5];
            this->state[3] ^= this->state[4];
            this->state[7] ^= this->state[3];
            this->state[5] ^= this->state[1];
            this->state[2] ^= this->state[0];
        } else {
            const UIntType state3_xor_state1 = rotr<w>(this->state[3], b);
            this->state[0] ^= state3_xor_state1;
            this->state[2] ^= this->state[0];
            this->state[1] ^= this->state[2] ^ this->state[0];
            this->state[1] = unshift_left_xor<w, a>(this->state[1]);
            UIntType t = this->state[1] << a;
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                t &= mask_result;
            }
            this->state[2] ^= t;
            this->state[3] = state3_xor_state1 ^ this->state[1];
        }

        return *this;
    }

    using common_engine<UIntType, w, xoshiro_engine>::operator++;
    using common_engine<UIntType, w, xoshiro_engine>::operator--;
    using common_engine<UIntType, w, xoshiro_engine>::operator==;

    auto operator==(const xoshiro_engine &rhs) const -> bool {
        return this->state == rhs.state;
    }

  private:
    std::array<UIntType, n> state;
};

using xoshiro512plus = xoshiro_engine<std::uint_fast64_t, 64, 8, 11, 21>;
using xoshiro512plusplus =
    xoshiro_engine<std::uint_fast64_t, 64, 8, 11, 21, 17>;
using xoshiro512starstar =
    xoshiro_engine<std::uint_fast64_t, 64, 8, 11, 21, 7, 5, 9>;

using xoshiro256plus =
    xoshiro_engine<std::uint_fast64_t, 64, 4, 17, 45, 0, 0, 0, 3, 2>;
using xoshiro256plusplus =
    xoshiro_engine<std::uint_fast64_t, 64, 4, 17, 45, 23, 0, 0, 3, 0>;
using xoshiro256starstar =
    xoshiro_engine<std::uint_fast64_t, 64, 4, 17, 45, 7, 5, 9>;

using xoshiro128plus =
    xoshiro_engine<std::uint_fast32_t, 32, 4, 9, 11, 0, 0, 0, 3, 2>;
using xoshiro128plusplus =
    xoshiro_engine<std::uint_fast32_t, 32, 4, 9, 11, 7, 0, 0, 3, 0>;
using xoshiro128starstar =
    xoshiro_engine<std::uint_fast32_t, 32, 4, 9, 11, 7, 5, 9>;

} // namespace adhoc

#endif // BRNG_XOSHIRO_ENGINE
