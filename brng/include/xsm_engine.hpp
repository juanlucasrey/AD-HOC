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

#ifndef BRNG_XSM_ENGINE
#define BRNG_XSM_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t r1, std::size_t r2,
          std::size_t r3, UIntType mult1, UIntType mult2 = mult1>
class xsm_engine final
    : public common_engine<UIntType, w,
                           xsm_engine<UIntType, w, r1, r2, r3, mult1, mult2>> {
    static constexpr UIntType mask_result = xsm_engine::max();

  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = w;

    static_assert(w <= std::numeric_limits<UIntType>::digits);

    static constexpr UIntType default_seed_1 = 0xf1ea5eed;
    static constexpr UIntType default_seed_2 = 0U;

    xsm_engine() : xsm_engine(default_seed_1) {}

    explicit xsm_engine(UIntType seed1, UIntType seed2 = default_seed_2) {

        if constexpr (w == 32) {
            this->lcg_adder_low = seed1 | 1U;
            this->lcg_adder_high = seed2;
            // the one bit of seed that did
            // not effect which cycle it's on
            this->lcg_high = this->lcg_adder_high + seed1 << (w - 1U);
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->lcg_high &= mask_result;
            }
        } else {
            this->lcg_adder_low = (seed1 << 1U) | 1U;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->lcg_adder_low &= mask_result;
            }
            // every bit of seed except the
            // highest bit gets used in the adder
            this->lcg_adder_high = (seed1 >> (w - 1U)) | (seed2 << 1U);
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->lcg_adder_high &= mask_result;
            }
            // and the highest bit of seed is used to
            // determine which end of the cycle we start at
            this->lcg_high =
                this->lcg_adder_high ^ ((seed2 >> (w - 1U)) << (w - 1U));
        }
        this->lcg_low = this->lcg_adder_low;

        this->operator++();
        this->operator++();
    }

    inline auto operator*() const -> value_type {
        UIntType result = this->lcg_high_prev + this->lcg_low_prev;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            result &= mask_result;
        }
        result = this->lcg_high_prev ^ rotl<w>(result, r1);
        UIntType tmp = result + this->lcg_adder_high;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            tmp &= mask_result;
        }
        result ^= rotl<w>(tmp, r2);
        result *= mult1;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            result &= mask_result;
        }

        tmp = result + this->lcg_high;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            tmp &= mask_result;
        }
        result ^= rotl<w>(tmp, r3);
        result *= mult2;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            result &= mask_result;
        }
        result ^= result >> r3;
        return result;
    };

    inline auto operator++() -> xsm_engine & {
        this->lcg_low_prev = this->lcg_low;
        this->lcg_high_prev = this->lcg_high;

        UIntType const tmp_lcg = this->lcg_low + this->lcg_adder_high;
        this->lcg_low += this->lcg_adder_low;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->lcg_low &= mask_result;
        }
        this->lcg_high +=
            tmp_lcg + ((this->lcg_low < this->lcg_adder_low) ? 1 : 0);
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->lcg_high &= mask_result;
        }

        return *this;
    }

    inline auto operator--() -> xsm_engine & {
        this->lcg_low = this->lcg_low_prev;
        this->lcg_high = this->lcg_high_prev;

        bool carry = this->lcg_low_prev < this->lcg_adder_low;
        this->lcg_low_prev -= this->lcg_adder_low;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->lcg_low_prev &= mask_result;
        }
        this->lcg_high_prev -=
            this->lcg_low_prev + this->lcg_adder_high + (carry ? 1 : 0);
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->lcg_high_prev &= mask_result;
        }

        return *this;
    }

    using common_engine<UIntType, w, xsm_engine>::operator++;
    using common_engine<UIntType, w, xsm_engine>::operator--;

    auto operator==(const xsm_engine &rhs) const -> bool {
        return (this->lcg_low_prev == rhs.lcg_low_prev) &&
               (this->lcg_high_prev == rhs.lcg_high_prev) &&
               (this->lcg_low == rhs.lcg_low) &&
               (this->lcg_high == rhs.lcg_high) &&
               (this->lcg_adder_low == rhs.lcg_adder_low) &&
               (this->lcg_adder_high == rhs.lcg_adder_high);
    }

  private:
    UIntType lcg_low_prev{0};
    UIntType lcg_high_prev{0};
    UIntType lcg_low{0};
    UIntType lcg_high{0};
    UIntType lcg_adder_low{0};
    UIntType lcg_adder_high{0};
};

using xsm32 =
    xsm_engine<std::uint_fast32_t, 32, 9, 19, 16, 0xD251CF2D, 0x299529B5>;
using xsm64 =
    xsm_engine<std::uint_fast64_t, 64, 16, 40, 32, 0xA3EC647659359ACDULL>;

} // namespace adhoc

#endif // BRNG_XSM_ENGINE
