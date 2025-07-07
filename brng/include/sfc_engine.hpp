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

#ifndef BRNG_SFC_ENGINE
#define BRNG_SFC_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/unshift.hpp"

#include <cstddef>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t p, std::size_t q,
          std::size_t r>
class sfc_engine final
    : public common_engine<UIntType, w, sfc_engine<UIntType, w, p, q, r>> {
  public:
    using value_type = UIntType;

    static_assert(p < w);
    static_assert(q < w);
    static_assert(r < w);

    static constexpr UIntType mask_result = sfc_engine::max();
    static constexpr UIntType default_seed = 0xcafef00dbeef5eedULL;

    sfc_engine() : sfc_engine(default_seed) {}

    explicit sfc_engine(UIntType seed) : sfc_engine(seed, seed, seed) {}

    explicit sfc_engine(UIntType seed1, UIntType seed2, UIntType seed3 = 0,
                        UIntType seed4 = 1, std::size_t iters = 12)
        : a_(seed1 & mask_result), b_(seed2 & mask_result),
          c_(seed3 & mask_result), d_(seed4 & mask_result),
          result_cache(a_ + b_ + d_) {

        for (std::size_t i = 0; i < iters; ++i) {
            this->operator++();
        }
    }

    inline auto operator*() const -> value_type { return this->result_cache; }

    inline auto operator++() -> sfc_engine & {
        ++this->d_;
        this->a_ = this->b_ ^ (this->b_ >> q);
        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            this->a_ &= mask_result;
        }

        this->b_ = this->c_ + (this->c_ << r);
        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            this->b_ &= mask_result;
        }
        this->c_ = rotl<w>(this->c_, p) + this->result_cache;
        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            this->c_ &= mask_result;
        }
        this->result_cache = this->a_ + this->b_ + this->d_;
        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            this->result_cache &= mask_result;
        }
        return *this;
    }

    inline auto operator--() -> sfc_engine & {
        UIntType c_prev = unshift_left_plus<w, r>(this->b_);
        this->b_ = unshift_right_xor<w, q>(this->a_);
        this->result_cache = this->c_ - rotl<w>(c_prev, p);
        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            this->result_cache &= mask_result;
        }
        this->c_ = c_prev;
        this->a_ = this->result_cache - this->b_ - --this->d_;
        if constexpr (w != std::numeric_limits<UIntType>::digits) {
            this->a_ &= mask_result;
        }
        return *this;
    }

    using common_engine<UIntType, w, sfc_engine>::operator++;
    using common_engine<UIntType, w, sfc_engine>::operator--;
    using common_engine<UIntType, w, sfc_engine>::operator==;

    auto operator==(const sfc_engine &rhs) const -> bool {
        return (this->a_ == rhs.a_) && (this->b_ == rhs.b_) &&
               (this->c_ == rhs.c_) && (this->d_ == rhs.d_);
    }

  private:
    UIntType a_, b_, c_, d_;
    UIntType result_cache;
};

using sfc64a = sfc_engine<std::uint_fast64_t, 64, 24, 11, 3>;
using sfc64b = sfc_engine<std::uint_fast64_t, 64, 25, 12, 3>;
using sfc64 = sfc64a;

using sfc32a = sfc_engine<std::uint_fast32_t, 32, 21, 9, 3>;
using sfc32b = sfc_engine<std::uint_fast32_t, 32, 15, 8, 3>;
using sfc32c = sfc_engine<std::uint_fast32_t, 32, 25, 8, 3>;
using sfc32 = sfc32a;

using sfc16a = sfc_engine<std::uint_fast16_t, 16, 4, 3, 2>;
using sfc16b = sfc_engine<std::uint_fast16_t, 16, 6, 5, 2>;
using sfc16c = sfc_engine<std::uint_fast16_t, 16, 4, 5, 3>;
using sfc16d = sfc_engine<std::uint_fast16_t, 16, 6, 5, 3>;
using sfc16e = sfc_engine<std::uint_fast16_t, 16, 7, 5, 3>;
using sfc16f = sfc_engine<std::uint_fast16_t, 16, 7, 3, 2>;
using sfc16 = sfc16d;

using sfc8 = sfc_engine<std::uint_fast8_t, 8, 3, 2, 1>;

} // namespace adhoc

#endif // BRNG_SFC_ENGINE
