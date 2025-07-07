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

#ifndef BRNG_WELL19937a
#define BRNG_WELL19937a

#include "tools/circular_buffer.hpp"
#include "tools/mask.hpp"
#include "tools/unshift.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, bool Tempering = false> class WELL19937a final {

    // this short step ensures consistency when the full_cycle is completed,
    // making simul.forward[full_cycle - 2] = simul.backward[1]
    void init_consistent() {
        this->operator()<true>();
        this->operator()<false>();
    }

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = 32U;
    static constexpr UIntType default_seed = 12345U;
    static_assert(word_size <= std::numeric_limits<UIntType>::digits);

    WELL19937a() : WELL19937a(default_seed) {}

    explicit WELL19937a(result_type value) {

        // values cannot all be 0
        if (value == 0) {
            value = default_seed;
        }

        this->state.data().fill(value);
        this->init_consistent();
    }

    template <class SeedSeq> explicit WELL19937a(SeedSeq &seq) {
        std::array<std::uint32_t, decltype(this->state)::size()>
            generated_sequence;
        seq.generate(generated_sequence.begin(), generated_sequence.end());

        std::transform(generated_sequence.begin(), generated_sequence.end(),
                       this->state.data().begin(), [](std::uint32_t v) {
                           return static_cast<UIntType>(v);
                       });

        // values cannot all be 0
        if (std::all_of(this->state.data().begin(), this->state.data().end(),
                        [](UIntType i) { return i == 0; })) {
            this->state.data().fill(default_seed);
        }
        this->init_consistent();
    }

    template <bool FwdDirection = true> inline auto operator()() -> UIntType {
        constexpr auto lower_mask = mask<UIntType>(31);
        constexpr auto upper_mask = ~lower_mask;
        constexpr auto global_mask = WELL19937a::max();

        UIntType result;
        if constexpr (FwdDirection) {
            auto const VRm1 = this->state.template at<623>();
            auto const VRm2 = this->state.template at<622>();
            auto const VM0 = this->state.at();
            auto const VM1 = this->state.template at<70>();
            auto const VM2 = this->state.template at<179>();
            auto const VM3 = this->state.template at<449>();
            auto &newV1 = this->state.at();
            auto &newV0 = this->state.template at<623>();

            this->cache[0] = (VRm1 & upper_mask) | (VRm2 & lower_mask);
            this->cache[1] = (VM0 ^ (VM0 << 25)) ^ (VM1 ^ (VM1 >> 27));
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                this->cache[1] &= global_mask;
            }
            UIntType const z2 = (VM2 >> 9) ^ (VM3 ^ (VM3 >> 1));
            newV1 = this->cache[1] ^ z2;
            newV0 = this->cache[0] ^ (this->cache[1] ^ (this->cache[1] << 9)) ^
                    (z2 ^ (z2 << 21)) ^ (newV1 ^ (newV1 >> 21));
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                newV0 &= global_mask;
            }

            --this->state;
            result = this->state.at();
        } else {
            result = this->state.at();
            ++this->state;

            auto &VM0 = this->state.at();
            auto const VM1 = this->state.template at<70>();

            UIntType const VM0Scrambled = this->cache[1] ^ (VM1 ^ (VM1 >> 27));
            VM0 = unshift_left_xor<word_size, 25U>(VM0Scrambled);
            UIntType const cache0_prev = this->cache[0];

            auto const newV1P = this->state.template at<1>();
            auto const VM2P = this->state.template at<180>();
            auto const VM3P = this->state.template at<450>();
            UIntType const z2P = (VM2P >> 9) ^ (VM3P ^ (VM3P >> 1));

            this->cache[1] = newV1P ^ z2P;
            this->cache[0] = VM0 ^ (this->cache[1] ^ (this->cache[1] << 9)) ^
                             (z2P ^ (z2P << 21)) ^ (newV1P ^ (newV1P >> 21));

            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                this->cache[0] &= global_mask;
            }

            auto &VRm1 = this->state.template at<623>();
            VRm1 = (cache0_prev & upper_mask) | (this->cache[0] & lower_mask);
        }

        if constexpr (Tempering) {
            constexpr UIntType TEMPERB = 0xE46E1700U;
            constexpr UIntType TEMPERC = 0x9B868000U;
            result = result ^ ((result << 7) & TEMPERB);
            result = result ^ ((result << 15) & TEMPERC);
        }
        return result;
    }

    static constexpr auto min() -> UIntType {
        return static_cast<UIntType>(0U);
    }

    static constexpr auto max() -> UIntType {
        return mask<UIntType>(word_size);
    }

    auto operator==(const WELL19937a &rhs) const -> bool {
        if (this == &rhs) {
            return true;
        }
        return this->state == rhs.state && this->cache == rhs.cache;
    }

    auto operator!=(const WELL19937a &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    circular_buffer<UIntType, 624> state{};
    std::array<UIntType, 2> cache{};
};

} // namespace adhoc

#endif // BRNG_WELL19937a
