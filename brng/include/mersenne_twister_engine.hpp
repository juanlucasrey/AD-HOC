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

#ifndef BRNG_MERSENNE_TWISTER_ENGINE
#define BRNG_MERSENNE_TWISTER_ENGINE

#include "tools/circular_buffer.hpp"
#include "tools/mask.hpp"
#include "tools/seed_seq_filler.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t n, std::size_t m,
          std::size_t r, UIntType a, std::size_t u, UIntType d, std::size_t s,
          UIntType b, std::size_t t, UIntType c, std::size_t l, UIntType f>
class mersenne_twister_engine final {
    static_assert(1 <= m && m <= n);
    static_assert(w >= 3);
    static_assert(w >= r);
    static_assert(w >= u);
    static_assert(w >= s);
    static_assert(w >= t);
    static_assert(w >= l);
    static_assert(w <= std::numeric_limits<UIntType>::digits);

    // this short step ensures consistency when the full_cycle is completed,
    // making simul.forward[full_cycle - 2] = simul.backward[1]
    void init_consistent() {
        this->operator()<true>();
        this->operator()<false>();
    }

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t state_size = n;
    static constexpr std::size_t shift_size = m;
    static constexpr std::size_t mask_bits = r;
    static constexpr UIntType xor_mask = a;
    static constexpr std::size_t tempering_u = u;
    static constexpr UIntType tempering_d = d;
    static constexpr std::size_t tempering_s = s;
    static constexpr UIntType tempering_b = b;
    static constexpr std::size_t tempering_t = t;
    static constexpr UIntType tempering_c = c;
    static constexpr std::size_t tempering_l = l;
    static constexpr UIntType initialization_multiplier = f;
    static constexpr UIntType default_seed = 5489U;

    mersenne_twister_engine() : mersenne_twister_engine(default_seed) {}

    explicit mersenne_twister_engine(result_type value) {
        constexpr auto mask = mersenne_twister_engine::max();
        this->state.at() = value;
        ++this->state;
        for (std::size_t j = 1; j < n; j++, ++this->state) {
            auto const prev = this->state.template at<n - 1>();
            this->state.at() =
                (f * (prev ^ (prev >> (w - 2))) + static_cast<UIntType>(j)) &
                mask;
        }

        this->init_consistent();
    }

    template <class SeedSeq> explicit mersenne_twister_engine(SeedSeq &seq) {

        seed_seq_filler<w, n> seq_filler(seq);
        seq_filler.template generate<UIntType>(state.data().begin(),
                                               state.data().end());

        this->init_consistent();
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        UIntType result;

        constexpr auto lower_mask = mask<UIntType>(r);
        constexpr auto upper_mask = ~lower_mask;

        if constexpr (FwdDirection) {
            auto const next_val = this->state.template at<1>();
            this->cache =
                (this->state.at() & upper_mask) | (next_val & lower_mask);
            this->state.at() = this->state.template at<m>() ^
                               (this->cache >> 1) ^ ((next_val & 1U) * a);

            result = this->state.at();
            ++this->state;
        } else {
            --this->state;
            result = this->state.at();

            auto const prev_cache = this->cache;

            this->cache = this->state.template at<n - 1>() ^
                          this->state.template at<m - 1>();
            UIntType const tmp = this->cache >> (w - 1);
            this->cache = ((this->cache ^ (tmp * a)) << 1) | tmp;
            this->state.at() =
                (prev_cache & upper_mask) | (this->cache & lower_mask);
        }

        result ^= ((result >> u) & d);
        result ^= ((result << s) & b);
        result ^= ((result << t) & c);
        result ^= (result >> l);
        return result;
    }

    void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }

    static constexpr auto min() -> UIntType {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> UIntType { return mask<UIntType>(w); }

    auto operator==(const mersenne_twister_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->cache == rhs.cache);
    }

    auto operator!=(const mersenne_twister_engine &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    circular_buffer<UIntType, n> state{};
    UIntType cache{0};
};

using mt19937 =
    mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31, 0x9908b0df,
                            11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18,
                            1812433253>;

using mt19937_64 =
    mersenne_twister_engine<std::uint_fast64_t, 64, 312, 156, 31,
                            0xb5026f5aa96619e9, 29, 0x5555555555555555, 17,
                            0x71d67fffeda60000, 37, 0xfff7eee000000000, 43,
                            6364136223846793005>;

} // namespace adhoc

#endif // BRNG_MERSENNE_TWISTER_ENGINE
