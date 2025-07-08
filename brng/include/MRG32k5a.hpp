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

#ifndef BRNG_MRG32k5a
#define BRNG_MRG32k5a

#include "tools/common_engine.hpp"
#include "tools/modular_arithmetic.hpp"

#include <array>
#include <cstdint>

namespace adhoc {

// The original implementation has logic that purposely avoids outputing 0.
// However this adds some branching that is a little convoluted when compared
// with the purely modular operator logic. I prefer using Original = false, and
// a wrapper that ensures 0 is never output in the left boundary of a
// uniform_real_distribution by other means.
template <class UIntType, bool Original = false>
class MRG32k5a final
    : public common_engine<UIntType, 32, MRG32k5a<UIntType, Original>> {
  public:
    using value_type = UIntType;

    MRG32k5a() : MRG32k5a(DefaultSeed) {}

    explicit MRG32k5a(UIntType value) {

        if (value == 0) {
            value = DefaultSeed;
        }
        // The seeds for Cg[0], Cg[1], Cg[2] must be integers in [0, m1 - 1] and
        // not all 0.

        // The seeds for Cg[3], Cg[4], Cg[5] must be integers in [0, m2 - 1] and
        // not all 0.

        this->state.fill(value);

        this->state[0] %= mod1;
        this->state[1] %= mod1;
        this->state[2] %= mod1;
        this->state[3] %= mod1;
        this->state[4] %= mod1;
        this->state[5] %= mod2;
        this->state[6] %= mod2;
        this->state[7] %= mod2;
        this->state[8] %= mod2;
        this->state[9] %= mod2;
        this->operator++();
    }

    template <class SeedSeq> explicit MRG32k5a(SeedSeq &seq) {
        seq.generate(this->state.begin(), this->state.end());
        this->state[0] %= mod1;
        this->state[1] %= mod1;
        this->state[2] %= mod1;
        this->state[3] %= mod1;
        this->state[4] %= mod1;
        this->state[5] %= mod2;
        this->state[6] %= mod2;
        this->state[7] %= mod2;
        this->state[8] %= mod2;
        this->state[9] %= mod2;

        if (this->state[0] == 0 && this->state[1] == 0 && this->state[2] == 0 &&
            this->state[3] == 0 && this->state[4] == 0) {
            this->state[0] = DefaultSeed;
            this->state[1] = DefaultSeed;
            this->state[2] = DefaultSeed;
            this->state[3] = DefaultSeed;
            this->state[4] = DefaultSeed;
        }

        if (this->state[5] == 0 && this->state[6] == 0 && this->state[7] == 0 &&
            this->state[8] == 0 && this->state[9] == 0) {
            this->state[5] = DefaultSeed;
            this->state[6] = DefaultSeed;
            this->state[7] = DefaultSeed;
            this->state[8] = DefaultSeed;
            this->state[9] = DefaultSeed;
        }
        this->operator++();
    }

    inline auto operator*() const -> value_type {
        if constexpr (Original) {
            return this->state[4] == this->state[9]
                       ? static_cast<UIntType>(mod1)
                       : static_cast<UIntType>(
                             (this->state[4] + (mod1 - this->state[9])) % mod1);

        } else {
            return static_cast<UIntType>(
                (this->state[4] + (mod1 - this->state[9])) % mod1);
        }
    }

    inline auto operator++() -> MRG32k5a & {
        auto const p1 =
            (((a12 * this->state[3] + a14 * this->state[1]) % mod1) +
             a15 * this->state[0]) %
            mod1;
        this->state[0] = this->state[1];
        this->state[1] = this->state[2];
        this->state[2] = this->state[3];
        this->state[3] = this->state[4];
        this->state[4] = p1;

        auto const p2 =
            (((a21 * this->state[9] + a23 * this->state[7]) % mod2) +
             a25 * this->state[5]) %
            mod2;
        this->state[5] = this->state[6];
        this->state[6] = this->state[7];
        this->state[7] = this->state[8];
        this->state[8] = this->state[9];
        this->state[9] = p2;
        return *this;
    }

    inline auto operator--() -> MRG32k5a & {
        constexpr auto a25_inv = modular_multiplicative_inverse(mod2, a25);
        constexpr auto a21_inv = mod2 - a21;
        constexpr auto a23_inv = mod2 - a23;

        auto const p2 = (((((a21_inv * this->state[8]) % mod2) +
                           (a23_inv * this->state[6]) + this->state[9]) %
                          mod2) *
                         a25_inv) %
                        mod2;
        this->state[9] = this->state[8];
        this->state[8] = this->state[7];
        this->state[7] = this->state[6];
        this->state[6] = this->state[5];
        this->state[5] = p2;

        constexpr auto a15_inv = modular_multiplicative_inverse(mod1, a15);
        constexpr auto a12_inv = mod1 - a12;
        constexpr auto a14_inv = mod1 - a14;
        auto const p1 = (((((a12_inv * this->state[2]) % mod1) +
                           (a14_inv * this->state[0]) + this->state[4]) %
                          mod1) *
                         a15_inv) %
                        mod1;
        this->state[4] = this->state[3];
        this->state[3] = this->state[2];
        this->state[2] = this->state[1];
        this->state[1] = this->state[0];
        this->state[0] = p1;
        return *this;
    }

    static constexpr auto min() -> UIntType {
        if constexpr (Original) {
            return static_cast<UIntType>(1U);
        } else {
            return static_cast<UIntType>(0U);
        }
    }

    static constexpr auto max() -> UIntType {
        if constexpr (Original) {
            return static_cast<UIntType>(mod1);
        } else {
            return static_cast<UIntType>(mod1 - 1);
        }
    }

    using common_engine<UIntType, 32, MRG32k5a>::operator++;
    using common_engine<UIntType, 32, MRG32k5a>::operator--;
    using common_engine<UIntType, 32, MRG32k5a>::operator==;

    auto operator==(const MRG32k5a &rhs) const -> bool {
        return (this->state == rhs.state);
    }

  private:
    static constexpr UIntType DefaultSeed = 12345U;
    static constexpr std::uint64_t a12 = 1154721U;
    static constexpr std::uint64_t a14 = 1739991U;
    static constexpr std::uint64_t a15 = 4293840528U;
    static constexpr std::uint64_t a21 = 1776413U;
    static constexpr std::uint64_t a23 = 865203U;
    static constexpr std::uint64_t a25 = 4293293275U;
    static constexpr std::uint64_t mod1 = 4294949027U;
    static constexpr std::uint64_t mod2 = 4294934327U;

    std::array<std::uint64_t, 10> state{};
};

} // namespace adhoc

#endif // BRNG_MRG32k5a
