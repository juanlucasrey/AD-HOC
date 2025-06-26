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

#ifndef BRNG_TOOLS_LINEAR_FEEDBACK_SHIFT_REGISTER
#define BRNG_TOOLS_LINEAR_FEEDBACK_SHIFT_REGISTER

#include "mask.hpp"
#include "unshift.hpp"

#include <algorithm>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t k, std::size_t q,
          std::size_t s>
class linear_feedback_shift_register final {
  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t low_bits = word_size - k;
    static constexpr std::size_t shift = k - s;
    static constexpr UIntType default_seed =
        w >= 64 ? 123456789123456789ULL : 987654321U;

    static_assert(k < word_size);
    static_assert(s < k);
    static_assert((s + q) < word_size);

    static_assert(word_size <= std::numeric_limits<UIntType>::digits);
    static_assert(std::is_unsigned_v<UIntType>);

    linear_feedback_shift_register()
        : linear_feedback_shift_register(default_seed) {}

    explicit linear_feedback_shift_register(result_type value) : state(value) {

        constexpr auto global_mask = linear_feedback_shift_register::max();

        constexpr UIntType lowval = static_cast<UIntType>(1) << low_bits;
        this->state = std::max(this->state, lowval) & global_mask;

        // this short step ensures consistency when the full_cycle is completed,
        // making simul.forward[full_cycle - 2] = simul.backward[1]
        this->operator()<true>();
        this->operator()<false>();
    };

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {

            if constexpr (w == std::numeric_limits<UIntType>::digits) {
                // in this case we can save an & operation
                UIntType const b = ((this->state << q) ^ this->state) >> shift;
                constexpr auto m = mask<UIntType>(shift, word_size - k);
                this->state = ((this->state & m) << s) ^ b;

            } else {
                constexpr auto global_mask =
                    linear_feedback_shift_register::max();
                UIntType const b =
                    (((this->state << q) & global_mask) ^ this->state) >> shift;
                constexpr auto m = mask<UIntType>(shift, word_size - k);
                this->state = ((this->state & m) << s) ^ b;
            }

            return this->state;
        } else {
            auto const result = this->state;
            constexpr std::size_t shift_inv = word_size - shift;
            constexpr std::size_t start = word_size - s - q;
            constexpr std::size_t high_scr = word_size - s;

            constexpr auto maskb = mask<UIntType>(shift_inv);
            UIntType const b = this->state & maskb;
            this->state = (this->state ^ b) >> s;

            unshift_left_xor<word_size, q, start>(this->state, b << shift);

            UIntType const bP = this->state & maskb;
            UIntType z_prev_middleP = (this->state ^ bP) >> s;

            unshift_left_xor<high_scr, q, start>(z_prev_middleP, bP << shift);
            UIntType const bPP =
                ((z_prev_middleP << q) ^ z_prev_middleP) >> shift;

            constexpr auto mask_low = mask<UIntType>(low_bits);
            this->state |= (mask_low & bPP);

            return result;
        }
    }

    static constexpr auto min() -> UIntType {
        return static_cast<UIntType>(0U);
    }

    static constexpr auto max() -> UIntType {
        return mask<UIntType>(word_size);
    }

    auto operator==(const linear_feedback_shift_register &rhs) const -> bool {
        return (this->state == rhs.state);
    }

    auto operator!=(const linear_feedback_shift_register &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    UIntType state{};
};

} // namespace adhoc

#endif // BRNG_TOOLS_LINEAR_FEEDBACK_SHIFT_REGISTER
