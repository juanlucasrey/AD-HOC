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

#ifndef BRNG_GJRAND_ENGINE
#define BRNG_GJRAND_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"

#include <cstddef>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, UIntType d_inc, std::size_t p,
          std::size_t q, std::size_t r>
class gjrand_engine final
    : public common_engine<UIntType, w,
                           gjrand_engine<UIntType, w, d_inc, p, q, r>> {
  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = w;

    static_assert(p < w);
    static_assert(q < w);
    static_assert(r < w);

    static constexpr UIntType mask_result = gjrand_engine::max();
    static constexpr UIntType default_seed_1 = 0xcafef00dbeef5eedULL;
    static constexpr UIntType default_seed_2 = 0ULL;
    static constexpr UIntType default_seed_3 = w <= 8    ? 201U
                                               : w <= 16 ? 2001U
                                               : w <= 32 ? 2000001U
                                                         : 2000001UL;

    gjrand_engine()
        : gjrand_engine(default_seed_1, default_seed_2, default_seed_3) {}

    explicit gjrand_engine(UIntType seed1, UIntType seed2, UIntType seed3)
        : a_(seed1 & mask_result), b_(seed2 & mask_result),
          c_(seed3 & mask_result), d_(0) {

        this->discard(15);
    }

    inline auto operator*() const -> value_type { return this->a_; }

    inline auto operator++() -> gjrand_engine & {
        this->b_ += this->c_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b_ &= mask_result;
        }
        this->a_ = rotl<w>(this->a_, p);
        this->c_ ^= this->b_;
        this->d_ += d_inc;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->d_ &= mask_result;
        }
        this->a_ += this->b_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a_ &= mask_result;
        }
        this->c_ = rotl<w>(this->c_, q);
        this->b_ ^= this->a_;
        this->a_ += this->c_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a_ &= mask_result;
        }
        this->b_ = rotl<w>(this->b_, r);
        this->c_ += this->a_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->c_ &= mask_result;
        }
        this->b_ += this->d_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b_ &= mask_result;
        }
        return *this;
    }

    inline auto operator--() -> gjrand_engine & {
        this->b_ -= this->d_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b_ &= mask_result;
        }

        this->c_ -= this->a_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->c_ &= mask_result;
        }
        this->b_ = rotr<w>(this->b_, r);
        this->a_ -= this->c_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a_ &= mask_result;
        }
        this->b_ ^= this->a_;
        this->c_ = rotr<w>(this->c_, q);
        this->a_ -= this->b_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a_ &= mask_result;
        }
        this->d_ -= d_inc;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->d_ &= mask_result;
        }
        this->c_ ^= this->b_;
        this->a_ = rotr<w>(this->a_, p);
        this->b_ -= this->c_;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b_ &= mask_result;
        }
        return *this;
    }

    using common_engine<UIntType, w, gjrand_engine>::operator++;
    using common_engine<UIntType, w, gjrand_engine>::operator--;
    using common_engine<UIntType, w, gjrand_engine>::operator==;

    auto operator==(const gjrand_engine &rhs) const -> bool {
        return (this->a_ == rhs.a_) && (this->b_ == rhs.b_) &&
               (this->c_ == rhs.c_) && (this->d_ == rhs.d_);
    }

  private:
    UIntType a_, b_, c_, d_;
};

using gjrand64 = gjrand_engine<std::uint64_t, 64U, 0x55aa96a5UL, 32U, 23U, 19U>;
using gjrand32 = gjrand_engine<std::uint32_t, 32U, 0x96a5U, 16U, 11U, 19U>;
using gjrand16 = gjrand_engine<std::uint16_t, 16U, 0x96a5U, 8U, 5U, 10U>;
using gjrand8 = gjrand_engine<std::uint8_t, 8U, 0x35U, 4U, 2U, 5U>;

} // namespace adhoc

#endif // BRNG_GJRAND_ENGINE
