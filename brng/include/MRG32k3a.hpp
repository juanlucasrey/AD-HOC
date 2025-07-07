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

#ifndef BRNG_MRG32k3a
#define BRNG_MRG32k3a

#include "tools/modular_arithmetic.hpp"

#include <array>
#include <cstdint>

namespace adhoc {

// The original implementation has logic that purposely avoids outputing 0.
// However this adds some branching that is a little convoluted when compared
// with the purely modular operator logic. I prefer using Original = false, and
// a wrapper that ensures 0 is never output in the left boundary of a
// uniform_real_distribution by other means.
template <class UIntType, bool Original = false> class MRG32k3a final {
    static_assert(std::numeric_limits<UIntType>::digits >= 32);
    static_assert(std::is_unsigned_v<UIntType>);

  public:
    using result_type = UIntType;

    MRG32k3a() : MRG32k3a(DefaultSeed) {}

    explicit MRG32k3a(result_type value) {

        if (value == 0) {
            value = DefaultSeed;
        }
        // The seeds for Cg[0], Cg[1], Cg[2] must be integers in [0, mod1 - 1]
        // and not all 0.

        // The seeds for Cg[3], Cg[4], Cg[5] must be integers in [0, mod2 - 1]
        // and not all 0.

        this->State.fill(value);

        this->State[0] %= mod1;
        this->State[1] %= mod1;
        this->State[2] %= mod1;
        this->State[3] %= mod2;
        this->State[4] %= mod2;
        this->State[5] %= mod2;
    }

    template <class SeedSeq> explicit MRG32k3a(SeedSeq &seq) {
        seq.generate(this->State.begin(), this->State.end());
        this->State[0] %= mod1;
        this->State[1] %= mod1;
        this->State[2] %= mod1;
        this->State[3] %= mod2;
        this->State[4] %= mod2;
        this->State[5] %= mod2;

        if (this->State[0] == 0 && this->State[1] == 0 && this->State[2] == 0) {
            this->State[0] = DefaultSeed;
            this->State[1] = DefaultSeed;
            this->State[2] = DefaultSeed;
        }

        if (this->State[3] == 0 && this->State[4] == 0 && this->State[5] == 0) {
            this->State[3] = DefaultSeed;
            this->State[4] = DefaultSeed;
            this->State[5] = DefaultSeed;
        }
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            auto const p1 =
                (((a13 * this->State[0]) % mod1) + (a12 * this->State[1])) %
                mod1;
            this->State[0] = this->State[1];
            this->State[1] = this->State[2];
            this->State[2] = p1;

            auto const p2 =
                (((a23 * this->State[3]) % mod2) + (a21 * this->State[5])) %
                mod2;
            this->State[3] = this->State[4];
            this->State[4] = this->State[5];
            this->State[5] = p2;

            if constexpr (Original) {
                return this->State[2] == this->State[5]
                           ? static_cast<UIntType>(mod1)
                           : static_cast<UIntType>(
                                 (this->State[2] + (mod1 - this->State[5])) %
                                 mod1);
            } else {
                return static_cast<UIntType>(
                    (this->State[2] + (mod1 - this->State[5])) % mod1);
            }
        } else {
            UIntType result = 0;
            if constexpr (Original) {
                result = this->State[2] == this->State[5]
                             ? static_cast<UIntType>(mod1)
                             : static_cast<UIntType>(
                                   (this->State[2] + (mod1 - this->State[5])) %
                                   mod1);
            } else {
                result = static_cast<UIntType>(
                    (this->State[2] + (mod1 - this->State[5])) % mod1);
            }

            constexpr auto a23_inv = modular_multiplicative_inverse(mod2, a23);
            constexpr auto a21_inv = mod2 - a21;

            auto const p2 =
                (((this->State[5] + a21_inv * this->State[4]) % mod2) *
                 a23_inv) %
                mod2;
            this->State[5] = this->State[4];
            this->State[4] = this->State[3];
            this->State[3] = p2;

            constexpr auto a13_inv = modular_multiplicative_inverse(mod1, a13);
            constexpr auto a12_inv = mod1 - a12;
            auto const p1 =
                (((this->State[2] + a12_inv * this->State[0]) % mod1) *
                 a13_inv) %
                mod1;
            this->State[2] = this->State[1];
            this->State[1] = this->State[0];
            this->State[0] = p1;

            return result;
        }
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

    auto operator==(const MRG32k3a &rhs) const -> bool {
        return (this->State == rhs.State);
    }

    auto operator!=(const MRG32k3a &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    static constexpr UIntType DefaultSeed = 12345U;
    static constexpr std::uint64_t a12 = 1403580U;
    static constexpr std::uint64_t a13 = 4294156359U;
    static constexpr std::uint64_t a21 = 527612U;
    static constexpr std::uint64_t a23 = 4293573854U;
    static constexpr std::uint64_t mod1 = 4294967087U;
    static constexpr std::uint64_t mod2 = 4294944443U;

    std::array<std::uint64_t, 6> State{};
};

} // namespace adhoc

#endif // BRNG_MRG32k3a
