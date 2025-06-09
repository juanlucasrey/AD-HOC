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

#include "tools/mask.hpp"

#include <cstddef>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, UIntType d_inc, std::size_t p,
          std::size_t q, std::size_t r>
class gjrand_engine final {
  public:
    using result_type = UIntType;

    static constexpr auto min() -> result_type {
        return static_cast<result_type>(0U);
    };

    static constexpr auto max() -> result_type { return mask<UIntType, w>(); };

    static_assert(w <= std::numeric_limits<UIntType>::digits);
    static_assert(p < w);
    static_assert(q < w);
    static_assert(r < w);

    static constexpr UIntType m = gjrand_engine::max();
    static constexpr UIntType default_seed_1 = 0xcafef00dbeef5eedULL;
    static constexpr UIntType default_seed_2 = 0ULL;
    static constexpr UIntType default_seed_3 = w <= 8    ? 201U
                                               : w <= 16 ? 2001U
                                               : w <= 32 ? 2000001U
                                                         : 2000001UL;

    gjrand_engine()
        : gjrand_engine(default_seed_1, default_seed_2, default_seed_3) {}

    explicit gjrand_engine(result_type seed1, result_type seed2,
                           result_type seed3)
        : a_(seed1 & m), b_(seed2 & m), c_(seed3 & m), d_(0) {

        for (std::size_t j = 0; j < 14; ++j) {
            this->operator()<true>();
        }
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        auto rotate = []<std::size_t k>(UIntType &x) {
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                x = ((x << k) & m) | (x >> (w - k));
            } else {
                x = (x << k) | (x >> (w - k));
            }
        };

        if constexpr (FwdDirection) {
            this->b_ += this->c_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->b_ &= m;
            }
            rotate.template operator()<p>(this->a_);
            this->c_ ^= this->b_;
            this->d_ += d_inc;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->d_ &= m;
            }
            this->a_ += this->b_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a_ &= m;
            }
            rotate.template operator()<q>(this->c_);
            this->b_ ^= this->a_;
            this->a_ += this->c_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a_ &= m;
            }
            rotate.template operator()<r>(this->b_);
            this->c_ += this->a_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c_ &= m;
            }
            this->b_ += this->d_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->b_ &= m;
            }

            return this->a_;
        } else {
            auto const result = this->a_;
            this->b_ -= this->d_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->b_ &= m;
            }

            this->c_ -= this->a_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c_ &= m;
            }
            rotate.template operator()<w - r>(this->b_);
            this->a_ -= this->c_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a_ &= m;
            }
            this->b_ ^= this->a_;
            rotate.template operator()<w - q>(this->c_);
            this->a_ -= this->b_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a_ &= m;
            }
            this->d_ -= d_inc;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->d_ &= m;
            }
            this->c_ ^= this->b_;
            rotate.template operator()<w - p>(this->a_);
            this->b_ -= this->c_;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->b_ &= m;
            }
            return result;
        }
    }

    template <bool FwdDirection = true> void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }

    auto operator==(const gjrand_engine &rhs) -> bool {
        return (this->a_ == rhs.a_) && (this->b_ == rhs.b_) &&
               (this->c_ == rhs.c_) && (this->d_ == rhs.d_);
    }

    auto operator!=(const gjrand_engine &rhs) -> bool {
        return !(this->operator==(rhs));
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
