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

#ifndef BRNG_WELL44497a
#define BRNG_WELL44497a

#include "tools/circular_buffer.hpp"
#include "tools/mask.hpp"
#include "tools/unshift.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, bool Tempering = false> class WELL44497a final {

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

    WELL44497a() : WELL44497a(default_seed) {}

    explicit WELL44497a(result_type value) {

        // values cannot all be 0
        if (value == 0) {
            value = default_seed;
        }

        this->state.data().fill(value);
        this->init_consistent();
    };

    template <class SeedSeq> explicit WELL44497a(SeedSeq &seq) {
        std::array<std::uint32_t, this->state.size()> generated_sequence;
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
        constexpr auto lower_mask = mask<UIntType, 15>();
        constexpr auto upper_mask = ~lower_mask;
        constexpr auto global_mask = WELL44497a::max();

        UIntType result;
        if constexpr (FwdDirection) {
            auto const VRm1 = this->state.template at<1390>();
            auto const VRm2 = this->state.template at<1389>();
            auto const VM0 = this->state.at();
            auto const VM1 = this->state.template at<23>();
            auto const VM2 = this->state.template at<481>();
            auto const VM3 = this->state.template at<229>();
            auto &newV1 = this->state.at();
            auto &newV0 = this->state.template at<1390>();

            this->cache[0] = (VRm1 & upper_mask) | (VRm2 & lower_mask);
            this->cache[1] = (VM0 ^ (VM0 << 24)) ^ (VM1 ^ (VM1 >> 30));
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                this->cache[1] &= global_mask;
            }
            UIntType z2 = (VM2 ^ (VM2 << 10)) ^ (VM3 << 26);
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                z2 &= global_mask;
            }
            newV1 = this->cache[1] ^ z2;
            UIntType MAT5 =
                ((z2 & 0x00020000U)
                     ? ((((z2 << 9) ^ (z2 >> (32 - 9))) & 0xfbffffffU) ^
                        0xb729fcecU)
                     : (((z2 << 9) ^ (z2 >> (32 - 9))) & 0xfbffffffU));
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                MAT5 &= global_mask;
            }

            newV0 = this->cache[0] ^ (this->cache[1] ^ (this->cache[1] >> 20)) ^
                    MAT5 ^ newV1;

            --this->state;
            result = this->state.at();
        } else {
            result = this->state.at();
            ++this->state;

            auto &VM0 = this->state.at();
            auto const VM1 = this->state.template at<23>();

            UIntType const VM0Scrambled = this->cache[1] ^ (VM1 ^ (VM1 >> 30));
            VM0 = unshift_left_xor<word_size, 24U>(VM0Scrambled);
            UIntType const cache0_prev = this->cache[0];

            auto const newV1P = this->state.template at<1>();
            auto const VM2P = this->state.template at<482>();
            auto const VM3P = this->state.template at<230>();
            UIntType z2P = (VM2P ^ (VM2P << 10)) ^ (VM3P << 26);
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                z2P &= global_mask;
            }
            UIntType MAT5P =
                ((z2P & 0x00020000U)
                     ? ((((z2P << 9) ^ (z2P >> (32 - 9))) & 0xfbffffffU) ^
                        0xb729fcecU)
                     : (((z2P << 9) ^ (z2P >> (32 - 9))) & 0xfbffffffU));
            if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
                MAT5P &= global_mask;
            }

            this->cache[1] = newV1P ^ z2P;
            this->cache[0] = VM0 ^ (this->cache[1] ^ (this->cache[1] >> 20)) ^
                             MAT5P ^ newV1P;

            auto &VRm1 = this->state.template at<1390>();
            VRm1 = (cache0_prev & upper_mask) | (this->cache[0] & lower_mask);
        }

        if constexpr (Tempering) {
            constexpr UIntType TEMPERB = 0x93dd1400U;
            constexpr UIntType TEMPERC = 0xfa118000U;
            result = result ^ ((result << 7) & TEMPERB);
            result = result ^ ((result << 15) & TEMPERC);
        }
        return result;
    }

    static constexpr auto min() -> UIntType {
        return static_cast<UIntType>(0U);
    }

    static constexpr auto max() -> UIntType {
        return mask<UIntType, word_size>();
    }

    auto operator==(const WELL44497a &rhs) const -> bool {
        if (this == &rhs) {
            return true;
        }
        return this->state == rhs.state && this->cache == rhs.cache;
    }

  private:
    circular_buffer<UIntType, 1391> state{};
    std::array<UIntType, 2> cache{};
};

} // namespace adhoc

#endif // BRNG_WELL44497a
