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

#ifndef BRNG_LINEAR_CONGRUENTIAL_ENGINE
#define BRNG_LINEAR_CONGRUENTIAL_ENGINE

#include "tools/mask.hpp"
#include "tools/modular_arithmetic.hpp"
#include "tools/uint128.hpp"

#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, UIntType a, UIntType c, UIntType m>
class linear_congruential_engine final {
    static_assert(m == 0 || a < m);
    static_assert(m == 0 || c < m);
    static_assert(m == 0 ? gcd(static_cast<UIntType>(2U), a) == 1U
                         : gcd(m, a) == 1U);
    static_assert(w <= std::numeric_limits<UIntType>::digits);

  public:
    using result_type = UIntType;

    static constexpr UIntType multiplier = a;
    static constexpr UIntType multiplier_inverse =
        m == 0 ? modular_multiplicative_inverse_pow2(w, a)
               : modular_multiplicative_inverse(m, a);
    static constexpr UIntType increment = c;
    static constexpr UIntType modulus = m;
    static constexpr UIntType default_seed = 1U;

    explicit linear_congruential_engine(result_type value = default_seed)
        : state(value) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        if constexpr (m == 0) {
            constexpr auto mask = linear_congruential_engine::max();
            if constexpr (FwdDirection) {
                if constexpr (w == std::numeric_limits<UIntType>::digits) {
                    this->state = a * this->state + c;
                } else {
                    this->state = (a * this->state + c) & mask;
                }
                return this->state;
            } else {
                const auto result = this->state;
                if constexpr (w == std::numeric_limits<UIntType>::digits) {
                    this->state = multiplier_inverse * (this->state - c);
                } else {
                    this->state =
                        (multiplier_inverse * (this->state - c)) & mask;
                }
                return result;
            }
        } else {
            using upgraded_type = std::conditional_t<
                w <= 8, std::uint_fast16_t,
                std::conditional_t<
                    w <= 16, std::uint_fast32_t,
                    std::conditional_t<w <= 32, std::uint_fast64_t, uint128>>>;

            constexpr auto m_64 = static_cast<upgraded_type>(m);
            if constexpr (c == 0) {
                if constexpr (FwdDirection) {
                    constexpr auto a_64 = static_cast<upgraded_type>(a);

                    this->state = static_cast<UIntType>(
                        (a_64 * static_cast<upgraded_type>(this->state)) %
                        m_64);
                    return this->state;
                } else {
                    constexpr auto a_inv_64 =
                        static_cast<upgraded_type>(multiplier_inverse);

                    const auto result = this->state;
                    this->state = static_cast<UIntType>(
                        (a_inv_64 * static_cast<upgraded_type>(this->state)) %
                        m_64);
                    return result;
                }
            } else {
                constexpr auto c_64 = static_cast<upgraded_type>(c);

                if constexpr (FwdDirection) {
                    constexpr auto a_64 = static_cast<upgraded_type>(a);

                    this->state = static_cast<UIntType>(
                        (a_64 * static_cast<upgraded_type>(this->state) +
                         c_64) %
                        m_64);
                    return this->state;
                } else {
                    constexpr auto a_inv_64 =
                        static_cast<upgraded_type>(multiplier_inverse);

                    constexpr auto c_inv_64 = static_cast<upgraded_type>(m - c);

                    const auto result = this->state;
                    this->state = static_cast<UIntType>(
                        (a_inv_64 *
                         (static_cast<upgraded_type>(this->state) + c_inv_64)) %
                        m_64);
                    return result;
                }
            }
        }
    }

    void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }

    static constexpr auto min() -> UIntType { return c == 0U ? 1U : 0U; }
    static constexpr auto max() -> UIntType {
        return m == 0 ? mask<UIntType>(w) : m - 1U;
    }

    auto operator==(const linear_congruential_engine &rhs) const -> bool {
        return (this->state == rhs.state);
    }

  private:
    UIntType state;
};

using minstd_rand0 =
    linear_congruential_engine<std::uint_fast32_t, 32, 16807, 0, 2147483647>;

using minstd_rand =
    linear_congruential_engine<std::uint_fast32_t, 32, 48271, 0, 2147483647>;

using ZX81 = linear_congruential_engine<std::uint_fast16_t, 16, 75, 74, 0>;

using ranqd1 =
    linear_congruential_engine<std::uint_fast32_t, 32, 1664525, 1013904223, 0>;

using RANDU =
    linear_congruential_engine<std::uint_fast32_t, 32, 65539, 0, 2147483648>;

using Borland =
    linear_congruential_engine<std::uint_fast32_t, 32, 22695477, 1, 2147483648>;

using Newlib =
    linear_congruential_engine<std::uint_fast64_t, 64, 6364136223846793005, 1,
                               9223372036854775808U>;

} // namespace adhoc

#endif // BRNG_LINEAR_CONGRUENTIAL_ENGINE
