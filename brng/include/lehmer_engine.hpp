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

#ifndef BRNG_LEHMER_ENGINE
#define BRNG_LEHMER_ENGINE

#include "tools/common_engine.hpp"
#include "tools/mask.hpp"
#include "tools/modular_arithmetic.hpp"
#include "tools/uint128.hpp"

#include <cstddef>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, class UIntType_int, std::size_t w_int,
          UIntType_int multiplier>
class lehmer_engine final
    : public common_engine<
          UIntType, w,
          lehmer_engine<UIntType, w, UIntType_int, w_int, multiplier>> {
  public:
    using value_type = UIntType;

    static_assert(w_int <= std::numeric_limits<UIntType_int>::digits);
    static_assert(w <= w_int);

    static constexpr UIntType mask_result = lehmer_engine::max();
    static constexpr UIntType_int m_int = mask<UIntType_int>(w_int);
    static constexpr UIntType default_seed = 0x9f57c403d06c42fcUL;

    lehmer_engine() : lehmer_engine(default_seed) {}

    explicit lehmer_engine(UIntType seed) : state_((seed & mask_result) | 1) {
        this->operator++();
    }

    inline auto operator*() const -> value_type {
        if constexpr (w_int == w) {
            return this->state_;
        } else {
            return static_cast<UIntType>(this->state_ >> (w_int - w));
        }
    }

    inline auto operator++() -> lehmer_engine & {
        this->state_ *= multiplier;

        if constexpr (w_int != std::numeric_limits<UIntType_int>::digits) {
            this->state_ &= m_int;
        }
        return *this;
    }

    inline auto operator--() -> lehmer_engine & {
        static constexpr UIntType_int multiplier_inverse =
            modular_multiplicative_inverse_pow2(w_int, multiplier);

        this->state_ *= multiplier_inverse;
        if constexpr (w_int != std::numeric_limits<UIntType_int>::digits) {
            this->state_ &= m_int;
        }
        return *this;
    }

    using common_engine<UIntType, w, lehmer_engine>::operator++;
    using common_engine<UIntType, w, lehmer_engine>::operator--;
    using common_engine<UIntType, w, lehmer_engine>::operator==;

    auto operator==(const lehmer_engine &rhs) const -> bool {
        return (this->state_ == rhs.state_);
    }

  private:
    UIntType_int state_;
};

using mcg128 = lehmer_engine<std::uint64_t, 64, uint128, 128,
                             (uint128(5017888479014934897ULL) << 64) +
                                 2747143273072462557ULL>;

using mcg128_fast =
    lehmer_engine<std::uint64_t, 64, uint128, 128, 0xda942042e4dd58b5ULL>;

} // namespace adhoc

#endif // BRNG_LEHMER_ENGINE
