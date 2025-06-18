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

#include "tools/mask.hpp"

#include <algorithm>
#include <array>
#include <cstdint>

namespace adhoc {

template <class UIntType> class arc4_engine final {
    static constexpr unsigned int BITS = 8;
    static constexpr unsigned int BOXES = 1U << BITS;

    void addrandom(const UIntType *dat, size_t datlen) {
        UIntType si;

        --this->index1;
        for (size_t n = 0; n < BOXES; ++n) {
            this->index1 = (this->index1 + 1) % BOXES;
            si = this->state[this->index1];
            this->index2 = (this->index2 + si + dat[n % datlen]) % BOXES;
            this->state[this->index1] = this->state[this->index2];
            this->state[this->index2] = si;
        }
        this->index2 = this->index1;

        // Throw away the first N bytes of output, as suggested in the
        // paper "Weaknesses in the Key Scheduling Algorithm of RC4"
        // by Fluher, Mantin, and Shamir.  N=1024 is based on
        // suggestions in the paper "(Not So) Random Shuffles of RC4"
        // by Ilya Mironov.
        for (size_t n = 0; n < 1024; ++n) {
            this->getval();
        }
    }

  public:
    using result_type = UIntType;
    static constexpr std::size_t word_size = 32;

    arc4_engine() : arc4_engine(123) {}

    explicit arc4_engine(result_type value) {};

    template <class SeedSeq> explicit arc4_engine(SeedSeq &seq) {
        for (size_t n = 0; n < BOXES; ++n) {
            this->state[n] = n;
        }

        constexpr size_t MAX_POSSIBLE_BOXES = 1U << BITS;
        constexpr bool IS_EFFICIENT = MAX_POSSIBLE_BOXES == BOXES;
        constexpr size_t BOXES_PER_INT32 =
            IS_EFFICIENT ? (sizeof(uint32_t) * CHAR_BIT) / BITS : 1;
        constexpr size_t SEED_BOXES = BOXES / 2;
        constexpr size_t SCRIBBLE_SIZE =
            (SEED_BOXES + BOXES_PER_INT32 - 1) / BOXES_PER_INT32;

        std::array<std::uint32_t, SCRIBBLE_SIZE> scribble;
        seq.generate(scribble.begin(), scribble.end());

        UIntType seeds[SEED_BOXES];
        for (size_t i = 0; i < SEED_BOXES; ++i) {
            seeds[i] = (scribble[i / BOXES_PER_INT32] >>
                        (BITS * (i % BOXES_PER_INT32))) %
                       BOXES;
        }
        this->addrandom(seeds, SEED_BOXES);
    }

    template <bool FwdDirection = true> inline auto getval() -> UIntType {
        if constexpr (FwdDirection) {
            this->index1 = (this->index1 + 1) % BOXES;
            this->index2 = (this->index2 + this->state[this->index1]) % BOXES;
            std::swap(this->state[this->index1], this->state[this->index2]);
            return this->state[(this->state[this->index2] +
                                this->state[this->index1]) %
                               BOXES];
        } else {
            auto result = this->state[(this->state[this->index2] +
                                       this->state[this->index1]) %
                                      BOXES];
            std::swap(this->state[this->index2], this->state[this->index1]);

            this->index2 =
                (this->index2 + (BOXES - this->state[this->index1])) % BOXES;
            this->index1 = (this->index1 + (BOXES - 1)) % BOXES;

            return result;
        }
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        result_type val = 0;
        constexpr unsigned int RESULT_BITS = sizeof(result_type) * CHAR_BIT;

        if constexpr (FwdDirection) {
            for (unsigned int shift = 0; shift < RESULT_BITS; shift += BITS) {
                val |= result_type(getval()) << shift;
            }
        } else {
            for (unsigned int shift = RESULT_BITS; shift;) {
                shift -= BITS;
                val |= result_type(getval<false>()) << shift;
            }
        }

        return val;
    }

    static constexpr auto min() -> UIntType {
        return static_cast<UIntType>(0U);
    }

    static constexpr auto max() -> UIntType {
        return mask<UIntType>(word_size);
    }

    auto operator==(const arc4_engine &rhs) const -> bool {
        return (this->index1 == rhs.index1) && (this->index2 == rhs.index2) &&
               (this->state == rhs.state);
    }

    auto operator!=(const arc4_engine &rhs) -> bool {
        return !(this->operator==(rhs));
    }

  private:
    UIntType index1 = 0;
    UIntType index2 = 0;
    std::array<UIntType, BOXES> state{};
};

using arc4_rand32 = arc4_engine<std::uint32_t>;
using arc4_rand64 = arc4_engine<std::uint64_t>;
using arc4_rand = arc4_rand32;

} // namespace adhoc

#endif // BRNG_ARC4_ENGINE
