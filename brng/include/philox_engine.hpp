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

#ifndef BRNG_PHILOX_ENGINE
#define BRNG_PHILOX_ENGINE

#include "tools/mask.hpp"
#include "tools/uint128.hpp"

#include <array>
#include <cstdint>
#include <limits>

namespace adhoc {

template <typename UIntType, std::size_t w, std::size_t n, std::size_t r,
          UIntType... consts>
class philox_engine {
  private:
    static_assert(sizeof...(consts) == n);
    static_assert(n == 2 || n == 4);
    static_assert(0 < r);
    static_assert(0 < w && w <= std::numeric_limits<UIntType>::digits);

    inline void increase_counter() {
        constexpr auto mask = philox_engine::max();
        std::size_t i = 0;
        do {
            this->X[i] = (this->X[i] + 1) & mask;
            ++i;
        } while (i < n && !this->X[i - 1]);
    }

    inline void decrease_counter() {
        constexpr auto mask = philox_engine::max();
        std::size_t i = 0;
        do {
            this->X[i] = (this->X[i] - 1) & mask;
            ++i;
        } while (i < n && (this->X[i - 1] == mask));
    }

    // this short step ensures consistency when the full_cycle is completed,
    // making simul.forward[full_cycle - 2] = simul.backward[1]
    void init_consistent() {
        this->decrease_counter();
        this->generate();
        this->increase_counter();
    }

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t word_count = n;
    static constexpr std::size_t round_count = r;
    static constexpr result_type default_seed = 20111115U;

    static constexpr auto min() -> result_type {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> result_type { return mask<UIntType, w>(); }

    philox_engine() : philox_engine(default_seed) {}

    explicit philox_engine(result_type value) {
        constexpr auto mask = philox_engine::max();
        this->K[0] = value & mask;
        this->init_consistent();
    }

    template <class SeedSeq> explicit philox_engine(SeedSeq &seq) {
        constexpr std::size_t p = ((w - 1) / 32) + 1;
        constexpr std::size_t n_half = n / 2;

        std::array<std::uint_least32_t, n_half * p> a;
        seq.generate(a.begin(), a.end());

        constexpr auto mask = philox_engine::max();
        auto iter = a.begin();
        for (std::size_t n_idx = 0; n_idx < n_half; n_idx++) {
            UIntType val = 0;
            for (std::size_t p_idx = 0; p_idx < p; ++p_idx) {
                val += static_cast<UIntType>(*iter++) << 32 * p_idx;
            }
            this->K[n_idx] = val & mask;
        }
        this->init_consistent();
    }

    void set_counter(const std::array<result_type, n> &c) {
        constexpr auto mask = philox_engine::max();
        for (std::size_t i = 0; i < n; ++i) {
            this->X[i] = c[i] & mask;
        }
        this->j = n - 1;
        this->decrease_counter();
        this->generate();
        this->increase_counter();
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            ++this->j;
            if (this->j == n) {
                this->generate();
                this->increase_counter();
                this->j = 0;
            }
            return Y[this->j];
        } else {
            const result_type result = Y[this->j];
            if (this->j == 0) {
                this->j = n;
                this->decrease_counter();
                this->decrease_counter();
                this->generate();
                this->increase_counter();
            }
            --this->j;
            return result;
        }
    }

    inline void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()<true>();
        }
    }

    auto operator==(const philox_engine &rhs) const -> bool {
        return (this->X == rhs.X) && (this->K == rhs.K) && (this->Y == rhs.Y) &&
               (this->j == rhs.j);
    }

  private:
    template <std::unsigned_integral U>
    inline auto mulhilo(U a, U b) -> std::pair<U, U> {

        using upgraded_type = std::conditional_t<
            w <= 8, std::uint_fast16_t,
            std::conditional_t<
                w <= 16, std::uint_fast32_t,
                std::conditional_t<w <= 32, std::uint_fast64_t, uint128>>>;

        const upgraded_type ab =
            static_cast<upgraded_type>(a) * static_cast<upgraded_type>(b);
        return {static_cast<U>(ab >> w),
                static_cast<U>(ab) & philox_engine::max()};
    }

    inline void generate() {
        constexpr auto mask = philox_engine::max();
        constexpr std::array<UIntType, n> consts_arr{consts...};
        if constexpr (n == 2) {
            result_type S0 = this->X[0];
            result_type S1 = this->X[1];
            result_type K0 = this->K[0];
            for (size_t i = 0; i < round_count; ++i) {
                auto [hi, lo] = this->mulhilo(S0, consts_arr[0]);
                S0 = hi ^ K0 ^ S1;
                S1 = lo;
                K0 = (K0 + consts_arr[1]) & mask;
            }
            this->Y[0] = S0;
            this->Y[1] = S1;

        } else if constexpr (n == 4) {
            result_type S0 = this->X[0];
            result_type S1 = this->X[1];
            result_type S2 = this->X[2];
            result_type S3 = this->X[3];
            result_type K0 = this->K[0];
            result_type K1 = this->K[1];
            for (size_t i = 0; i < r; ++i) {
                auto [hi0, lo0] = this->mulhilo(S2, consts_arr[0]);
                auto [hi2, lo2] = this->mulhilo(S0, consts_arr[2]);
                S0 = hi0 ^ K0 ^ S1;
                S1 = lo0;
                S2 = hi2 ^ K1 ^ S3;
                S3 = lo2;
                K0 = (K0 + consts_arr[1]) & mask;
                K1 = (K1 + consts_arr[3]) & mask;
            }
            this->Y[0] = S0;
            this->Y[1] = S1;
            this->Y[2] = S2;
            this->Y[3] = S3;
        }
    }

    std::array<UIntType, n> X{0};
    std::array<UIntType, n / 2> K{0};
    std::array<UIntType, n> Y{0};
    std::size_t j{n - 1};
};

using philox4x32 = philox_engine<std::uint_fast32_t, 32, 4, 10, 0xCD9E8D57,
                                 0x9E3779B9, 0xD2511F53, 0xBB67AE85>;

using philox4x64 =
    philox_engine<std::uint_fast64_t, 64, 4, 10, 0xCA5A826395121157,
                  0x9E3779B97F4A7C15, 0xD2E7470EE14C6C93, 0xBB67AE8584CAA73B>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_PHILOX_ENGINE
