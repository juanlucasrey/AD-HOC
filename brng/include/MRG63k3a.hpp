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

#ifndef BRNG_MRG63k3a
#define BRNG_MRG63k3a

#include "tools/common_engine.hpp"
#include "tools/modular_arithmetic.hpp"
#include "tools/seed_seq_filler.hpp"
#include "tools/uint128.hpp"

#include <array>
#include <cstdint>

namespace adhoc {

// The original implementation has logic that purposely avoids outputing 0.
// However this adds some branching that is a little convoluted when compared
// with the purely modular operator logic. I prefer using Original = false, and
// a wrapper that ensures 0 is never output in the left boundary of a
// uniform_real_distribution by other means.
template <class UIntType, bool Use128 = false, bool Original = false>
class MRG63k3a final
    : public common_engine<UIntType, 64, MRG63k3a<UIntType, Use128, Original>> {
  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = 64;

    MRG63k3a() : MRG63k3a(DefaultSeed) {}

    explicit MRG63k3a(UIntType value) {

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

    template <class SeedSeq> explicit MRG63k3a(SeedSeq &seq) {

        seed_seq_filler<word_size, 6> seq_filler(seq);
        seq_filler.template generate<UIntType>(this->state.begin(),
                                               this->state.end());
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
        if constexpr (Use128) {
            if constexpr (Original) {
                return this->state[2] == this->state[5]
                           ? static_cast<UIntType>(mod1)
                           : static_cast<UIntType>(
                                 (this->state[2] + (mod1 - this->state[5])) %
                                 mod1);
            } else {
                return static_cast<UIntType>(
                    (this->state[2] + (mod1 - this->state[5])) % mod1);
            }
        } else {
            if constexpr (Original) {
                return (this->state[2] > this->state[5])
                           ? this->state[2] - this->state[5]
                           : (this->state[2] + (mod1 - this->state[5]));
            } else {
                return (this->state[2] >= this->state[5])
                           ? this->state[2] - this->state[5]
                           : (this->state[2] + (mod1 - this->state[5]));
            }
        }
    }

    inline auto operator++() -> MRG63k3a & {
        if constexpr (Use128) {
            static constexpr uint128 a13 = mod1 - a13n;
            static_assert(a13 == 9223372033672665121);
            static constexpr uint128 a12l = a12;
            static constexpr uint128 m1l = mod1;

            auto const p1 = ((a12l * static_cast<uint128>(this->state[1])) +
                             (a13 * static_cast<uint128>(this->state[0]))) %
                            m1l;
            this->state[0] = this->state[1];
            this->state[1] = this->state[2];
            this->state[2] = static_cast<std::uint64_t>(p1);

            static constexpr uint128 a23 = mod2 - a23n;
            static constexpr uint128 a21l = a21;
            static constexpr uint128 m2l = mod2;
            auto const p2 = ((a21l * static_cast<uint128>(this->state[5])) +
                             (a23 * static_cast<uint128>(this->state[3]))) %
                            m2l;
            this->state[3] = this->state[4];
            this->state[4] = this->state[5];
            this->state[5] = static_cast<std::uint64_t>(p2);
        } else {
            std::uint64_t h = 0;
            std::uint64_t v1 = 0;
            std::uint64_t v2 = 0;
            std::uint64_t p13 = 0;
            std::uint64_t p1 = 0;
            std::uint64_t p23 = 0;
            std::uint64_t p2 = 0;
            constexpr auto q13 = mod1 / a13n;
            constexpr auto r13 = mod1 % a13n;
            h = this->state[0] / q13;
            v1 = a13n * (this->state[0] - h * q13);
            v2 = h * r13;
            p13 = v1 >= v2 ? v1 - v2 : v1 + (mod1 - v2);

            constexpr auto q12 = mod1 / a12;
            constexpr auto r12 = mod1 % a12;
            h = this->state[1] / q12;
            v1 = a12 * (this->state[1] - h * q12);
            v2 = h * r12;
            p1 = v1 >= v2 ? v1 - v2 : v1 + mod1 - v2;

            p1 = p1 >= p13 ? p1 - p13 : p1 + (mod1 - p13);

            this->state[0] = this->state[1];
            this->state[1] = this->state[2];
            this->state[2] = p1;

            constexpr auto q23 = mod2 / a23n;
            constexpr auto r23 = mod2 % a23n;
            h = this->state[3] / q23;
            v1 = a23n * (this->state[3] - h * q23);
            v2 = h * r23;
            p23 = v1 >= v2 ? v1 - v2 : v1 + (mod2 - v2);

            constexpr auto q21 = mod2 / a21;
            constexpr auto r21 = mod2 % a21;
            h = this->state[5] / q21;
            v1 = a21 * (this->state[5] - h * q21);
            v2 = h * r21;
            p2 = v1 >= v2 ? v1 - v2 : v1 + mod2 - v2;

            p2 = p2 >= p23 ? p2 - p23 : p2 + (mod2 - p23);

            this->state[3] = this->state[4];
            this->state[4] = this->state[5];
            this->state[5] = p2;
        }

        return *this;
    }

    inline auto operator--() -> MRG63k3a & {
        if constexpr (Use128) {
            static constexpr auto a23 = mod2 - a23n;
            static constexpr uint128 m2l = mod2;
            static constexpr uint128 a21l_inv = mod2 - a21;

            constexpr uint128 a23_inv =
                modular_multiplicative_inverse(mod2, a23);

            auto const p2 =
                (((static_cast<uint128>(this->state[5]) +
                   a21l_inv * static_cast<uint128>(this->state[4])) %
                  m2l) *
                 a23_inv) %
                m2l;
            this->state[5] = this->state[4];
            this->state[4] = this->state[3];
            this->state[3] = static_cast<std::uint64_t>(p2);

            static constexpr auto a13 = mod1 - a13n;
            static constexpr uint128 m1l = mod1;
            static constexpr uint128 a12l_inv = mod1 - a12;

            constexpr uint128 a13_inv =
                modular_multiplicative_inverse(mod1, a13);

            auto const p1 =
                (((static_cast<uint128>(this->state[2]) +
                   (a12l_inv * static_cast<uint128>(this->state[0]))) %
                  m1l) *
                 a13_inv) %
                m1l;
            this->state[2] = this->state[1];
            this->state[1] = this->state[0];
            this->state[0] = static_cast<std::uint64_t>(p1);
        } else {
            std::uint64_t h = 0;
            std::uint64_t v1 = 0;
            std::uint64_t v2 = 0;
            std::uint64_t p1 = 0;
            std::uint64_t p2 = 0;

            constexpr auto q21 = mod2 / a21;
            constexpr auto r21 = mod2 % a21;
            h = this->state[4] / q21;
            v1 = a21 * (this->state[4] - h * q21);
            v2 = h * r21;
            p2 = v1 >= v2 ? v1 - v2 : v1 + (mod2 - v2);

            p2 = this->state[5] >= p2 ? this->state[5] - p2
                                      : this->state[5] + (mod2 - p2);

            static constexpr auto a23 = mod2 - a23n;
            constexpr auto a23_inv = modular_multiplicative_inverse(mod2, a23);

            constexpr std::uint64_t a23_inv_1 = 2843350197U;
            constexpr std::uint64_t a23_inv_2 = 3103994195U;
            static_assert(a23_inv == a23_inv_1 * a23_inv_2);

            constexpr auto q3 = mod2 / a23_inv_1;
            constexpr auto r3 = mod2 % a23_inv_1;
            h = p2 / q3;
            v1 = a23_inv_1 * (p2 - h * q3);
            v2 = h * r3;
            p2 = v1 >= v2 ? v1 - v2 : v1 + (mod2 - v2);

            constexpr auto q4 = mod2 / a23_inv_2;
            constexpr auto r4 = mod2 % a23_inv_2;
            h = p2 / q4;
            v1 = a23_inv_2 * (p2 - h * q4);
            v2 = h * r4;
            p2 = v1 >= v2 ? v1 - v2 : v1 + (mod2 - v2);

            this->state[5] = this->state[4];
            this->state[4] = this->state[3];
            this->state[3] = p2;

            constexpr auto q12 = mod1 / a12;
            constexpr auto r12 = mod1 % a12;
            h = this->state[0] / q12;
            v1 = a12 * (this->state[0] - h * q12);
            v2 = h * r12;
            p1 = v1 >= v2 ? v1 - v2 : v1 + (mod1 - v2);

            p1 = this->state[2] >= p1 ? this->state[2] - p1
                                      : this->state[2] + (mod1 - p1);

            static constexpr auto a13 = mod1 - a13n;
            constexpr auto a13_inv = modular_multiplicative_inverse(mod1, a13);
            static_assert(a13_inv == 8929150787489966579);

            constexpr std::uint64_t a13_inv_1 = 2375717617U;
            constexpr std::uint64_t a13_inv_2 = 3758506787U;
            static_assert(a13_inv == a13_inv_1 * a13_inv_2);

            constexpr auto q5 = mod1 / a13_inv_1;
            constexpr auto r5 = mod1 % a13_inv_1;
            h = p1 / q5;
            v1 = a13_inv_1 * (p1 - h * q5);
            v2 = h * r5;
            p1 = v1 >= v2 ? v1 - v2 : v1 + (mod1 - v2);

            constexpr auto q6 = mod1 / a13_inv_2;
            constexpr auto r6 = mod1 % a13_inv_2;
            h = p1 / q6;
            v1 = a13_inv_2 * (p1 - h * q6);
            v2 = h * r6;
            p1 = v1 >= v2 ? v1 - v2 : v1 + (mod1 - v2);

            this->state[2] = this->state[1];
            this->state[1] = this->state[0];
            this->state[0] = p1;
        }
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

    using common_engine<UIntType, 64, MRG63k3a>::operator++;
    using common_engine<UIntType, 64, MRG63k3a>::operator--;
    using common_engine<UIntType, 64, MRG63k3a>::operator==;

    auto operator==(const MRG63k3a &rhs) const -> bool {
        return (this->state == rhs.state);
    }

  private:
    static constexpr UIntType DefaultSeed = 123456789U;
    static constexpr std::uint64_t a12 = 1754669720U;
    static constexpr std::uint64_t a13n = 3182104042U;
    static constexpr std::uint64_t a21 = 31387477935U;
    static constexpr std::uint64_t a23n = 6199136374U;
    static constexpr std::uint64_t mod1 = 9223372036854769163U;
    static constexpr std::uint64_t mod2 = 9223372036854754679U;

    std::array<std::uint64_t, 6> state{};
};

} // namespace adhoc

#endif // BRNG_MRG63k3a
