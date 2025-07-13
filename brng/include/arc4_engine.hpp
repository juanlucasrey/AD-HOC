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

#ifndef BRNG_ARC4_ENGINE
#define BRNG_ARC4_ENGINE

#include "tools/common_engine.hpp"

#include <algorithm>
#include <array>
#include <cstdint>

namespace adhoc {

template <class UIntType>
class arc4_engine final
    : public common_engine<UIntType, 8, arc4_engine<UIntType>> {
    void init() {
        // Throw away the first N bytes of output, as suggested in the
        // paper "Weaknesses in the Key Scheduling Algorithm of RC4"
        // by Fluher, Mantin, and Shamir.  N=1024 is based on
        // suggestions in the paper "(Not So) Random Shuffles of RC4"
        // by Ilya Mironov.
        this->discard(1025);
    }

  public:
    using value_type = UIntType;

    static constexpr std::size_t word_size = 8;
    static constexpr UIntType mask_result = arc4_engine::max();

    arc4_engine() : arc4_engine(123) {}

    explicit arc4_engine(UIntType value) {
        this->state[0] = value;
        this->init();
    };

    template <class SeedSeq> explicit arc4_engine(SeedSeq &seq) {
        constexpr unsigned int state_size = 1U << word_size;
        for (size_t n = 0; n < this->state.size(); ++n) {
            this->state[n] = n;
        }

        constexpr size_t BOXES_PER_INT32 =
            (sizeof(uint32_t) * CHAR_BIT) / word_size;
        constexpr size_t SEED_BOXES = state_size / 2;
        constexpr size_t SCRIBBLE_SIZE =
            (SEED_BOXES + BOXES_PER_INT32 - 1) / BOXES_PER_INT32;

        std::array<std::uint32_t, SCRIBBLE_SIZE> scribble;
        seq.generate(scribble.begin(), scribble.end());

        UIntType seeds[SEED_BOXES];
        for (size_t i = 0; i < SEED_BOXES; ++i) {
            seeds[i] = (scribble[i / BOXES_PER_INT32] >>
                        (word_size * (i % BOXES_PER_INT32))) %
                       this->state.size();
        }

        UIntType si;

        --this->index1;
        for (size_t n = 0; n < this->state.size(); ++n) {
            this->index1 = (this->index1 + 1) & mask_result;
            si = this->state[this->index1];
            this->index2 =
                (this->index2 + si + seeds[n % SEED_BOXES]) & mask_result;
            this->state[this->index1] = this->state[this->index2];
            this->state[this->index2] = si;
        }
        this->index2 = this->index1;
        this->init();
    }

    inline auto operator*() const -> value_type {
        return this
            ->state[(this->state[this->index2] + this->state[this->index1]) &
                    mask_result];
    }

    inline auto operator++() -> arc4_engine & {
        this->index1 = (this->index1 + 1) & mask_result;
        this->index2 = (this->index2 + this->state[this->index1]) & mask_result;
        std::swap(this->state[this->index1], this->state[this->index2]);
        return *this;
    }

    inline auto operator--() -> arc4_engine & {
        std::swap(this->state[this->index2], this->state[this->index1]);
        this->index2 =
            (this->index2 + ((mask_result + 1) - this->state[this->index1])) &
            mask_result;
        this->index1 = (this->index1 + mask_result) & mask_result;
        return *this;
    }

    using common_engine<UIntType, 8, arc4_engine>::operator++;
    using common_engine<UIntType, 8, arc4_engine>::operator--;
    using common_engine<UIntType, 8, arc4_engine>::operator==;

    auto operator==(const arc4_engine &rhs) const -> bool {
        return (this->index1 == rhs.index1) && (this->index2 == rhs.index2) &&
               (this->state == rhs.state);
    }

  private:
    UIntType index1 = 0;
    UIntType index2 = 0;
    std::array<UIntType, (1U << word_size)> state{};
};

using arc4 = arc4_engine<std::uint8_t>;

} // namespace adhoc

#endif // BRNG_ARC4_ENGINE
