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

#ifndef BRNG_MRG31k3p
#define BRNG_MRG31k3p

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
class MRG31k3p final
    : public common_engine<UIntType, 32, MRG31k3p<UIntType, Original>> {
  public:
    using value_type = UIntType;

    MRG31k3p() : MRG31k3p(DefaultSeed) {}

    explicit MRG31k3p(UIntType value) {

        if (value == 0) {
            value = DefaultSeed;
        }
        // The seeds for Cg[0], Cg[1], Cg[2] must be integers in [0, mod1 - 1]
        // and not all 0.

        // The seeds for Cg[3], Cg[4], Cg[5] must be integers in [0, mod2 - 1]
        // and not all 0.

        this->state.fill(value);

        this->state[0] %= mod1;
        this->state[1] %= mod1;
        this->state[2] %= mod1;
        this->state[3] %= mod2;
        this->state[4] %= mod2;
        this->state[5] %= mod2;
        this->operator++();
    }

    template <class SeedSeq> explicit MRG31k3p(SeedSeq &seq) {
        seq.generate(this->state.begin(), this->state.end());
        this->state[0] %= mod1;
        this->state[1] %= mod1;
        this->state[2] %= mod1;
        this->state[3] %= mod2;
        this->state[4] %= mod2;
        this->state[5] %= mod2;

        if (this->state[0] == 0 && this->state[1] == 0 && this->state[2] == 0) {
            this->state[0] = DefaultSeed;
            this->state[1] = DefaultSeed;
            this->state[2] = DefaultSeed;
        }

        if (this->state[3] == 0 && this->state[4] == 0 && this->state[5] == 0) {
            this->state[3] = DefaultSeed;
            this->state[4] = DefaultSeed;
            this->state[5] = DefaultSeed;
        }
        this->operator++();
    }

    inline auto operator*() const -> value_type {
        if constexpr (Original) {
            auto const p1_mod = this->state[2] == 0 ? mod1 : this->state[2];
            auto const p2_mod = this->state[5] == 0 ? mod2 : this->state[5];
            return p1_mod == p2_mod ? static_cast<UIntType>(mod1)
                                    : static_cast<UIntType>(
                                          (p1_mod + (mod1 - p2_mod)) % mod1);

        } else {
            return static_cast<UIntType>(
                (this->state[2] + (mod1 - this->state[5])) % mod1);
        }
    }

    inline auto operator++() -> MRG31k3p & {
        auto const p1 = (a13 * this->state[0] + a12 * this->state[1]) % mod1;
        this->state[0] = this->state[1];
        this->state[1] = this->state[2];
        this->state[2] = p1;

        auto const p2 = (a23 * this->state[3] + a21 * this->state[5]) % mod2;
        this->state[3] = this->state[4];
        this->state[4] = this->state[5];
        this->state[5] = p2;
        return *this;
    }

    inline auto operator--() -> MRG31k3p & {
        constexpr auto a23_inv = modular_multiplicative_inverse(mod2, a23);
        constexpr auto a21_inv = mod2 - a21;

        auto const p2 =
            (((this->state[5] + a21_inv * this->state[4]) % mod2) * a23_inv) %
            mod2;
        this->state[5] = this->state[4];
        this->state[4] = this->state[3];
        this->state[3] = p2;

        constexpr auto a13_inv = modular_multiplicative_inverse(mod1, a13);
        constexpr auto a12_inv = mod1 - a12;
        auto const p1 =
            (((this->state[2] + a12_inv * this->state[0]) % mod1) * a13_inv) %
            mod1;
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

    using common_engine<UIntType, 32, MRG31k3p>::operator++;
    using common_engine<UIntType, 32, MRG31k3p>::operator--;
    using common_engine<UIntType, 32, MRG31k3p>::operator==;

    auto operator==(const MRG31k3p &rhs) const -> bool {
        return (this->state == rhs.state);
    }

  private:
    static constexpr UIntType DefaultSeed = 12345U;
    static constexpr std::uint64_t a12 = (static_cast<std::uint64_t>(1) << 22U);
    static constexpr std::uint64_t a13 =
        (static_cast<std::uint64_t>(1) << 7U) + 1U;
    static constexpr std::uint64_t a21 =
        (static_cast<std::uint64_t>(1U) << 15U);
    static constexpr std::uint64_t a23 = a21 + 1U;
    static constexpr std::uint64_t mod1 = 2147483647U;
    static constexpr std::uint64_t mod2 = 2147462579U;

    std::array<std::uint64_t, 6> state{};
};

} // namespace adhoc

#endif // BRNG_MRG31k3p
