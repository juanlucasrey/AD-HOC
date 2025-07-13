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

#ifndef BRNG_RARNS_ENGINE
#define BRNG_RARNS_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/unshift.hpp"

#include <array>
#include <cstddef>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t rot1, std::size_t rot2,
          std::size_t rot3, std::size_t rs>
class rarns_engine final
    : public common_engine<UIntType, w,
                           rarns_engine<UIntType, w, rot1, rot2, rot3, rs>> {
  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = w;

    static constexpr UIntType mask_result = rarns_engine::max();
    static constexpr UIntType default_seed_1 = 0xf1ea5eed;
    static constexpr UIntType default_seed_2 = 0xcafe5eed00000001ULL;
    static constexpr std::size_t default_seed_3 = 20U;

    rarns_engine() : rarns_engine(default_seed_1) {}

    explicit rarns_engine(UIntType seed1, UIntType seed2 = default_seed_2,
                          UIntType seed3 = default_seed_3,
                          std::size_t iters = 5)
        : state({seed1, seed2, seed3}) {

        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            for (auto &s : this->state) {
                s &= mask_result;
            }
        }

        if constexpr (w == 32) {
            this->state[2] = ~(this->state[0] + this->state[1]);
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                this->state[2] &= mask_result;
            }
        }
        for (int i = 0; i < iters; i++) {
            this->operator++();
        }
        this->operator++();
    }

    inline auto operator*() const -> value_type {
        UIntType result = this->state_prev[0] + this->state_prev[2];
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            result &= mask_result;
        }
        result = rotl<w>(result, rot1);
        result += this->state[0];
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            result &= mask_result;
        }
        return result;
    }

    inline auto operator++() -> rarns_engine & {
        this->state_prev = this->state;
        UIntType const old = this->state[0] >> rs;
        this->state[2] ^= this->state[0];
        this->state[0] ^= this->state[1];
        this->state[1] ^= this->state[2];
        this->state[0] ^= old;
        this->state[2] = rotl<w>(this->state[2], rot2);
        this->state[0] = rotl<w>(this->state[0], rot3);
        return *this;
    }

    inline auto operator--() -> rarns_engine & {
        this->state = this->state_prev;
        this->state_prev[0] = rotr<w>(this->state_prev[0], rot3);
        this->state_prev[2] = rotr<w>(this->state_prev[2], rot2);
        this->state_prev[1] ^= this->state_prev[2];
        this->state_prev[0] ^= this->state_prev[1];
        this->state_prev[0] = unshift_right_xor<w, rs>(this->state_prev[0]);
        this->state_prev[2] ^= this->state_prev[0];
        return *this;
    }

    using common_engine<UIntType, w, rarns_engine>::operator++;
    using common_engine<UIntType, w, rarns_engine>::operator--;
    using common_engine<UIntType, w, rarns_engine>::operator==;

    auto operator==(const rarns_engine &rhs) const -> bool {
        return (this->state == rhs.state) &&
               (this->state_prev == rhs.state_prev);
    }

  private:
    std::array<UIntType, 3> state{};
    std::array<UIntType, 3> state_prev{};
};

template <class UIntType>
using rarns64 = rarns_engine<UIntType, 64, 21, 3, 5, 13>;
template <class UIntType>
using rarns32 = rarns_engine<UIntType, 32, 10, 3, 5, 13>;
template <class UIntType>
using rarns16 = rarns_engine<UIntType, 16, 5, 7, 8, 3>;

} // namespace adhoc

#endif // BRNG_RARNS_ENGINE
