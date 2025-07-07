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

#ifndef BRNG_SPLITMIX_ENGINE
#define BRNG_SPLITMIX_ENGINE

#include "tools/common_engine.hpp"
#include "tools/modular_arithmetic.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::uint64_t m1, std::uint64_t m2,
          std::size_t p, std::size_t q, std::size_t r, std::uint64_t m3,
          std::uint64_t m4, std::size_t s, std::size_t t, std::size_t u>
class splitmix_engine final
    : public common_engine<
          UIntType, w,
          splitmix_engine<UIntType, w, m1, m2, p, q, r, m3, m4, s, t, u>> {
    static_assert(w <= 64);

  public:
    using value_type = UIntType;

    static constexpr std::uint64_t default_seed1 = 0xbad0ff1ced15ea5eUL;
    static constexpr std::uint64_t default_seed2 = 0x9e3779b97f4a7c15UL;

    explicit splitmix_engine(std::uint64_t seed = default_seed1,
                             std::uint64_t gamma = default_seed2)
        : seed_(seed), gamma_(gamma | 1U), // gamma NEEDS to be odd so that it
                                           // is coprime with 2^w, otherwise the
          // modular multiplicative inverse would not exist!
          gamma_inverse_(modular_multiplicative_inverse_pow2(64, gamma_)) {}

    inline auto operator*() const -> value_type {
        std::uint64_t result = this->seed_;
        result ^= result >> s;
        result *= m3;
        result ^= result >> t;
        result *= m4;
        if constexpr (u) {
            result ^= result >> u;
        }

        if constexpr (w == 64) {
            return static_cast<value_type>(result);
        } else {
            return static_cast<value_type>(result >> (64 - w));
        }
    }

    inline auto operator++() -> splitmix_engine & {
        this->seed_ += this->gamma_;
        return *this;
    }

    inline auto operator--() -> splitmix_engine & {
        this->seed_ -= this->gamma_;
        return *this;
    }

    template <bool FwdDirection = true> void discard(unsigned long long z) {
        if constexpr (FwdDirection) {
            seed_ += z * gamma_;
        } else {
            seed_ -= z * gamma_;
        }
    }

    auto operator-(const splitmix_engine &rhs) -> unsigned long long {
        if (this->gamma_ == rhs.gamma_) {
            return (seed_ - rhs.seed_) * this->gamma_inverse_;
        } else {
            return std::numeric_limits<unsigned long long>::max();
        }
    }

    using common_engine<UIntType, w, splitmix_engine>::operator++;
    using common_engine<UIntType, w, splitmix_engine>::operator--;
    using common_engine<UIntType, w, splitmix_engine>::operator==;

    auto operator==(const splitmix_engine &rhs) const -> bool {
        return (this->seed_ == rhs.seed_) && (this->gamma_ == rhs.gamma_) &&
               (this->gamma_inverse_ == rhs.gamma_inverse_);
    }

    auto split() -> splitmix_engine {
        std::uint64_t new_seed = this->operator*();
        this->operator++();
        std::uint64_t x = this->seed_;
        x ^= x >> p;
        x *= m1;
        x ^= x >> q;
        x *= m2;
        x ^= x >> r;

        // gamma NEEDS to be odd so that it is coprime with 2^w, otherwise the
        // modular multiplicative inverse would not exist!
        x |= 1UL;

        unsigned int n = std::popcount(x ^ (x >> 1U));
        std::uint64_t new_gamma = (n < 24) ? x ^ 0xaaaaaaaaaaaaaaaaUL : x;

        this->operator++();
        return splitmix_engine{new_seed, new_gamma};
    }

  private:
    std::uint64_t seed_;
    std::uint64_t gamma_;
    std::uint64_t gamma_inverse_;
};

using splitmix64 =
    splitmix_engine<std::uint64_t, 64, 0xff51afd7ed558ccdUL,
                    0xc4ceb9fe1a85ec53UL, 33, 33, 33, 0xbf58476d1ce4e5b9UL,
                    0x94d049bb133111ebUL, 30, 27, 31>;

using splitmix32 =
    splitmix_engine<std::uint32_t, 32, 0xff51afd7ed558ccdUL,
                    0xc4ceb9fe1a85ec53UL, 33, 33, 33, 0x62a9d9ed799705f5UL,
                    0xcb24d0a5c88c35b3UL, 33, 28, 0>;

} // namespace adhoc

#endif // BRNG_SPLITMIX_ENGINE
