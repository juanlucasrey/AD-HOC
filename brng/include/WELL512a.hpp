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

#ifndef BRNG_WELL512a
#define BRNG_WELL512a

#include "tools/circular_buffer.hpp"
#include "tools/common_engine.hpp"
#include "tools/unshift.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType>
class WELL512a final : public common_engine<UIntType, 32, WELL512a<UIntType>> {
  public:
    using value_type = UIntType;

    static constexpr std::size_t word_size = 32U;
    static constexpr UIntType default_seed = 12345U;
    static_assert(word_size <= std::numeric_limits<UIntType>::digits);

    WELL512a() : WELL512a(default_seed) {}

    explicit WELL512a(UIntType value) {

        // values cannot all be 0
        if (value == 0) {
            value = default_seed;
        }

        this->state.data().fill(value);
        this->operator++();
    }

    template <class SeedSeq> explicit WELL512a(SeedSeq &seq) {
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
        this->operator++();
    }

    inline auto operator*() const -> value_type { return this->state.at(); }

    inline auto operator++() -> WELL512a & {
        constexpr auto global_mask = WELL512a::max();
        auto const &VM0 = this->state.at();
        auto const &VM1 = this->state.template at<13>();
        auto const &VM2 = this->state.template at<9>();
        auto const &VRm1 = this->state.template at<15>();
        auto &newV0 = this->state.template at<15>();
        auto &newV1 = this->state.at();

        UIntType z1 = (VM0 ^ (VM0 << 16U)) ^ (VM1 ^ (VM1 << 15U));
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            z1 &= global_mask;
        }
        UIntType const z2 = (VM2 ^ (VM2 >> 11U));
        newV1 = z1 ^ z2;
        newV0 = (VRm1 ^ (VRm1 << 2U)) ^ (z1 ^ (z1 << 18U)) ^ (z2 << 28U) ^
                (newV1 ^ ((newV1 << 5U) & 0xDA442D24U));
        if constexpr (word_size != std::numeric_limits<UIntType>::digits) {
            newV0 &= global_mask;
        }

        --this->state;
        return *this;
    }

    inline auto operator--() -> WELL512a & {
        constexpr auto global_mask = WELL512a::max();
        ++this->state;

        auto &VM0 = this->state.at();
        auto const &VM1 = this->state.template at<13>();
        auto const &VM2 = this->state.template at<9>();
        auto &VRm1 = this->state.template at<15>();
        auto const &newV0 = this->state.template at<15>();
        auto const &newV1 = this->state.at();

        UIntType const z2 = (VM2 ^ (VM2 >> 11U));
        UIntType const z1 = newV1 ^ z2;
        UIntType const VRm1Scrambled = newV0 ^ (z1 ^ (z1 << 18U)) ^
                                       (z2 << 28U) ^
                                       (newV1 ^ ((newV1 << 5U) & 0xDA442D24U));
        UIntType const VM0Scrambled = z1 ^ (VM1 ^ (VM1 << 15U));
        VRm1 = unshift_left_xor<word_size, 2U>(VRm1Scrambled);
        VM0 = unshift_left_xor<word_size, 16U>(VM0Scrambled);
        return *this;
    }

    using common_engine<UIntType, 32, WELL512a>::operator++;
    using common_engine<UIntType, 32, WELL512a>::operator--;
    using common_engine<UIntType, 32, WELL512a>::operator==;

    auto operator==(const WELL512a &rhs) const -> bool {
        return this->state == rhs.state;
    }

  private:
    circular_buffer<UIntType, 16> state;
};

} // namespace adhoc

#endif // BRNG_WELL512a
