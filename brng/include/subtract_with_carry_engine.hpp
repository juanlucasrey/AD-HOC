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

#ifndef BRNG_SUBTRACT_WITH_CARRY_ENGINE
#define BRNG_SUBTRACT_WITH_CARRY_ENGINE

#include "linear_congruential_engine.hpp"
#include "tools/circular_buffer.hpp"
#include "tools/mask.hpp"

#include <algorithm>
#include <array>
#include <limits>

namespace adhoc {

// G. Marsaglia, A. Zaman, A new class of random number generators, The Annals
// of Applied Probability 1 (3) (1991) 462–480. URL
// http://www.jstor.org/stable/2959748

// on page 11 of the article (page 472 of the publication) it says: "Finding
// conditions which ensure that a seed vector is periodic is pretty much an
// acedemic exercise, done out of curiosity"

// the corrected=true template parameter ensures that our RNG is strictly
// periodic.

// Marsaglia et al. pointed at this when he wrote on page 3 (page 464 of the
// publication) "since f has an inverse [...] the sequence is strictly periodic"
// corrected = true ensures that f has an inverse and is therefore strictly
// periodic.

template <class UIntType, std::size_t w, std::size_t s, std::size_t r>
class subtract_with_carry_engine final {
  private:
    static constexpr std::size_t k = (w + 31U) / 32U;

    void init(std::array<std::uint_least32_t, r * k> const &seeds) {
        constexpr auto mask = subtract_with_carry_engine::max();

        auto iter = seeds.begin();
        auto &data = this->state.data();
        for (auto &entry : data) {
            UIntType val = 0;
            for (std::size_t ki = 0; ki < k; ++ki) {
                val += static_cast<UIntType>(*iter++) << 32 * ki;
            }
            entry = val & mask;
        }

        this->carry = (this->state.template at<long_lag - 1>() == 0);

        for (std::size_t j = 0; j < long_lag; ++j) {
            this->operator()<true>();
        }
        for (std::size_t j = 0; j < long_lag; ++j) {
            this->operator()<false>();
        }
    }

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t short_lag = s;
    static constexpr std::size_t long_lag = r;
    static constexpr std::uint_least32_t default_seed = 19780503U;

    static constexpr result_type modulus = static_cast<result_type>(1U) << w;

    static_assert(s >= 1U);
    static_assert(s < r);
    static_assert(w >= 1U);
    static_assert(w <= std::numeric_limits<UIntType>::digits);

    subtract_with_carry_engine() : subtract_with_carry_engine(0U) {}

    explicit subtract_with_carry_engine(result_type value) {
        linear_congruential_engine<std::uint_least32_t, 32, 40014U, 0U,
                                   2147483563U>
            e(static_cast<std::uint_least32_t>(value == 0U ? default_seed
                                                           : value));
        std::array<std::uint_least32_t, r * k> seeds;
        std::generate(seeds.begin(), seeds.end(), e);
        this->init(seeds);
    }

    template <class SeedSeq> explicit subtract_with_carry_engine(SeedSeq &seq) {
        std::array<std::uint_least32_t, r * k> a;
        seq.generate(a.begin(), a.end());
        this->init(a);
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        constexpr auto mask = subtract_with_carry_engine::max();

        if constexpr (FwdDirection) {
            const result_type &xs =
                this->state.template at<long_lag - short_lag>();
            result_type &xr = this->state.at();

            const result_type new_carry = this->carry == 0 ? xs < xr : xs <= xr;
            xr = (xs - xr - this->carry) & mask;
            this->carry = new_carry;
            ++this->state;
            return xr;
        } else {
            --this->state;

            // usingn circular buffer is slower because of modulo ops
            constexpr bool use_circular_buffer = false;

            if constexpr (use_circular_buffer) {
                const result_type &xs =
                    this->state.template at<long_lag - short_lag>();
                result_type &xr = this->state.at();

                if (xs != xr) {
                    std::size_t count = 0;
                    do {
                        --this->state;
                        count++;
                        // this termination is safe:
                        // -if it never terminates, it means all states are
                        // equal but if all states are equal then either temp ==
                        // 0 or temp == modulus in which case we would never be
                        // here.
                    } while (this->state.template at<long_lag - short_lag>() ==
                             this->state.at());

                    this->carry =
                        this->state.template at<long_lag - short_lag>() <
                        this->state.at();

                    this->state += count;
                }
                const UIntType result = xr;
                xr = (xs - xr - this->carry) & mask;
                return result;
            } else {
                auto const &state_data = this->state.data();
                std::size_t xri_prev = this->state.index();
                std::size_t xsi_prev = (xri_prev < short_lag)
                                           ? (xri_prev + long_lag - short_lag)
                                           : (xri_prev - short_lag);

                const result_type &xs = state_data[xsi_prev];
                result_type &xr = this->state.at();

                if (xs != xr) {
                    do {
                        xri_prev =
                            (xri_prev == 0) ? (long_lag - 1) : (xri_prev - 1);
                        xsi_prev =
                            (xsi_prev == 0) ? (long_lag - 1) : (xsi_prev - 1);
                        // this termination is safe:
                        // -if it never terminates, it means all states are
                        // equal but if all states are equal then either temp ==
                        // 0 or temp == modulus in which case we would never be
                        // here.
                    } while (state_data[xsi_prev] == state_data[xri_prev]);

                    this->carry = state_data[xsi_prev] < state_data[xri_prev];
                }

                const UIntType result = xr;
                xr = (xs - xr - this->carry) & mask;
                return result;
            }
        }
    }

    void discard(unsigned long long z) {
        for (unsigned long long j = 0; j < z; ++j) {
            this->operator()();
        }
    }

    static constexpr auto min() -> UIntType {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> UIntType { return mask<UIntType, w>(); }

    auto operator==(const subtract_with_carry_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->carry == rhs.carry);
    }

  private:
    circular_buffer<UIntType, long_lag> state{};
    UIntType carry{0};
};

using ranlux24_base =
    subtract_with_carry_engine<std::uint_fast32_t, 24, 10, 24>;

using ranlux48_base = subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12>;

} // namespace adhoc

#endif // BRNG_SUBTRACT_WITH_CARRY_ENGINE
