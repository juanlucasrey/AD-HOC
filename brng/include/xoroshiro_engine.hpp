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

#ifndef BRNG_XOROSHIRO_ENGINE
#define BRNG_XOROSHIRO_ENGINE

#include "tools/bit.hpp"
#include "tools/circular_buffer.hpp"
#include "tools/common_engine.hpp"
#include "tools/seed_seq_filler.hpp"

#include <cstddef>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t bits, std::size_t a,
          std::size_t b, std::size_t c, UIntType mult1 = 0,
          std::size_t orot = 0, UIntType mult2 = 0, bool rotate_state = false>
class xoroshiro_engine final
    : public common_engine<UIntType, w,
                           xoroshiro_engine<UIntType, w, bits, a, b, c, mult1,
                                            orot, mult2, rotate_state>> {
  public:
    using value_type = UIntType;

    static_assert(w <= std::numeric_limits<UIntType>::digits);
    static_assert(a < w);
    static_assert(b < w);
    static_assert(c < w);

    static constexpr UIntType size = bits / w;
    static constexpr UIntType mask_result = xoroshiro_engine::max();
    static constexpr UIntType default_seed_1 = 0xc1f651c67c62c6e0;
    static constexpr UIntType default_seed_2 = 0x30d89576f866ac9f;

    xoroshiro_engine() : xoroshiro_engine(default_seed_1, default_seed_2) {}

    explicit xoroshiro_engine(UIntType seed1, UIntType seed2) {
        this->state.at() = seed1 & mask_result;
        this->state.template at<1>() =
            (seed1 || seed2) ? (seed2 & mask_result) : static_cast<UIntType>(1);

        if constexpr (rotate_state) {
            ++this->state;
        }
    }

    template <class SeedSeq> explicit xoroshiro_engine(SeedSeq &seq) {
        seed_seq_filler<w, size> seq_filler(seq);
        seq_filler.template generate<UIntType>(this->state.data().begin(),
                                               this->state.data().end());
        if constexpr (rotate_state) {
            ++this->state;
        }
    }

    inline auto operator*() const -> value_type {
        UIntType result = 0;

        if constexpr (mult1 && orot) {
            result = this->state.at() * mult1;
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                result &= mask_result;
            }
            result = rotl<w>(result, orot) * mult2;
        } else if constexpr (mult1) {
            result = this->state.at() * mult1;
        } else if constexpr (static_cast<bool>(orot)) {
            result = this->state.at() + this->state.template at<size - 1>();
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                result &= mask_result;
            }
            if constexpr (rotate_state) {
                result =
                    rotl<w>(result, orot) + this->state.template at<size - 1>();
            } else {
                result = rotl<w>(result, orot) + this->state.at();
            }
        } else {
            result = this->state.at() + this->state.template at<size - 1>();
        }

        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            result &= mask_result;
        }

        return result;
    }

    inline auto operator++() -> xoroshiro_engine & {
        auto &s0 = this->state.at();
        auto &s1 = this->state.template at<size - 1>();
        s1 ^= s0;
        if constexpr (rotate_state) {
            auto const s1_prev = s1;
            s1 = rotl<w>(s0, a) ^ s1 ^ (s1 << b);
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                s1 &= mask_result;
            }
            s0 = rotl<w>(s1_prev, c);
            ++this->state;
        } else {
            s0 = rotl<w>(s0, a) ^ s1 ^ (s1 << b);
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                s0 &= mask_result;
            }
            s1 = rotl<w>(s1, c);
        }
        return *this;
    }

    inline auto operator--() -> xoroshiro_engine & {
        if constexpr (rotate_state) {
            --this->state;
        }
        auto &s0 = this->state.at();
        auto &s1 = this->state.template at<size - 1>();
        if constexpr (rotate_state) {
            auto const s1_prev = rotr<w>(s0, c);
            s0 = s1 ^ s1_prev ^ (s1_prev << b);
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                s0 &= mask_result;
            }
            s0 = rotr<w>(s0, a);
            s1 = s1_prev;
        } else {
            s1 = rotr<w>(s1, c);
            s0 = s0 ^ s1 ^ (s1 << b);
            if constexpr (w != std::numeric_limits<UIntType>::digits) {
                s0 &= mask_result;
            }
            s0 = rotr<w>(s0, a);
        }
        s1 ^= s0;
        return *this;
    }

    using common_engine<UIntType, w, xoroshiro_engine>::operator++;
    using common_engine<UIntType, w, xoroshiro_engine>::operator--;
    using common_engine<UIntType, w, xoroshiro_engine>::operator==;

    auto operator==(const xoroshiro_engine &rhs) const -> bool {
        return (this->state == rhs.state);
    }

  private:
    circular_buffer<UIntType, size> state{};
};

using xoroshiro64star =
    xoroshiro_engine<std::uint_fast32_t, 32, 64, 26, 9, 13, 0x9E3779BB>;
using xoroshiro64starstar =
    xoroshiro_engine<std::uint_fast32_t, 32, 64, 26, 9, 13, 0x9E3779BB, 5, 5>;

using xoroshiro128plus =
    xoroshiro_engine<std::uint_fast64_t, 64, 128, 24, 16, 37>;
using xoroshiro128plusplus =
    xoroshiro_engine<std::uint_fast64_t, 64, 128, 49, 21, 28, 0, 17>;
using xoroshiro128starstar =
    xoroshiro_engine<std::uint_fast64_t, 64, 128, 24, 16, 37, 5, 7, 9>;

using xoroshiro1024plusplus =
    xoroshiro_engine<std::uint_fast64_t, 64, 1024, 25, 27, 36, 0, 23, 0, true>;
using xoroshiro1024star = xoroshiro_engine<std::uint_fast64_t, 64, 1024, 25, 27,
                                           36, 0x9e3779b97f4a7c13, 0, 0, true>;
using xoroshiro1024starstar =
    xoroshiro_engine<std::uint_fast64_t, 64, 1024, 25, 27, 36, 5, 7, 9, true>;

} // namespace adhoc

#endif // BRNG_XOROSHIRO_ENGINE
