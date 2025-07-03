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

#ifndef BRNG_THREEFRY_ENGINE
#define BRNG_THREEFRY_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t n, std::size_t r,
          UIntType key_add, int... rots>
class threefry_engine
    : public common_engine<
          UIntType, w, threefry_engine<UIntType, w, n, r, key_add, rots...>> {
  private:
    static_assert(sizeof...(rots) == n * 4);
    static_assert(n == 2 || n == 4);
    static_assert(0 < r);
    static constexpr UIntType mask_result = threefry_engine::max();

    inline void increase_counter() {
        std::size_t i = 0;
        do {
            this->X[i] = (this->X[i] + 1) & mask_result;
            ++i;
        } while (i < n && !this->X[i - 1]);
    }

    inline void decrease_counter() {
        std::size_t i = 0;
        do {
            this->X[i] = (this->X[i] - 1) & mask_result;
            ++i;
        } while (i < n && (this->X[i - 1] == mask_result));
    }

    void init() {
        this->K.back() = key_add;
        for (std::size_t i = 0; i < n; ++i) {
            this->K.back() ^= this->K[i];
        }

        this->operator++();
    }

    template <std::size_t I = 0> void generate_n2() {
        constexpr std::array<std::size_t, sizeof...(rots)> rots_a{rots...};
        constexpr std::size_t rot_idx = I % (n * 4);
        this->Y[0] += this->Y[1];
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->Y[0] &= mask_result;
        }
        this->Y[1] = rotl<w>(this->Y[1], rots_a[rot_idx]);
        this->Y[1] ^= this->Y[0];
        if constexpr ((I % 4) == (4 - 1)) {
            constexpr std::size_t add = (I / 4) + 1;
            constexpr std::size_t idxk1 = add % (n + 1);
            constexpr std::size_t idxk2 = (add + 1) % (n + 1);
            this->Y[0] += this->K[idxk1];
            this->Y[1] += this->K[idxk2];
            this->Y[1] += add;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->Y[0] &= mask_result;
                this->Y[1] &= mask_result;
            }
        }
        if constexpr ((I + 1) < r) {
            generate_n2<I + 1>();
        }
    }

    template <std::size_t I = 0> void generate_n4() {
        constexpr std::array<std::size_t, sizeof...(rots)> rots_a{rots...};
        constexpr std::size_t rot_idx1 = (I * 2) % (n * 4);
        constexpr std::size_t rot_idx2 = rot_idx1 + 1;
        constexpr std::size_t idx1 = ((I % 2) * 2) + 1;
        constexpr std::size_t idx2 = (((I + 1) % 2) * 2) + 1;
        this->Y[0] += this->Y[idx1];
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->Y[0] &= mask_result;
        }
        this->Y[idx1] = rotl<w>(this->Y[idx1], rots_a[rot_idx1]);
        this->Y[idx1] ^= this->Y[0];
        this->Y[2] += this->Y[idx2];
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->Y[2] &= mask_result;
        }
        this->Y[idx2] = rotl<w>(this->Y[idx2], rots_a[rot_idx2]);
        this->Y[idx2] ^= this->Y[2];
        if constexpr ((I % 4) == (4 - 1)) {
            constexpr std::size_t add = (I / 4) + 1;
            constexpr std::size_t idxk1 = add % (n + 1);
            constexpr std::size_t idxk2 = (add + 1) % (n + 1);
            constexpr std::size_t idxk3 = (add + 2) % (n + 1);
            constexpr std::size_t idxk4 = (add + 3) % (n + 1);
            this->Y[0] += this->K[idxk1];
            this->Y[1] += this->K[idxk2];
            this->Y[2] += this->K[idxk3];
            this->Y[3] += this->K[idxk4];
            this->Y[3] += add;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->Y[0] &= mask_result;
                this->Y[1] &= mask_result;
                this->Y[2] &= mask_result;
                this->Y[3] &= mask_result;
            }
        }
        if constexpr ((I + 1) < r) {
            generate_n4<I + 1>();
        }
    }

    inline void generate() {
        for (std::size_t i = 0; i < n; ++i) {
            this->Y[i] = this->X[i] + this->K[i];
        }

        if constexpr (n == 2) {
            this->generate_n2();
        } else {
            this->generate_n4();
        }
    }

  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = w;
    static constexpr std::size_t word_count = n;
    static constexpr std::size_t round_count = r;
    static constexpr UIntType default_seed = 20111115U;

    template <class... value_type>
    explicit threefry_engine(value_type... values) {
        static_assert(sizeof...(values) <= n);
        if constexpr (sizeof...(values)) {
            auto make_vector =
                []<class... Args>(Args... args) -> std::vector<int> {
                return {args...};
            };

            auto inputs = make_vector(values...);

            std::copy(inputs.begin(), inputs.end(), this->K.begin());

            for (auto &val : this->K) {
                val &= mask_result;
            }

        } else {
            this->K.front() = default_seed & mask_result;
        }

        this->init();
    }

    template <class SeedSeq> explicit threefry_engine(SeedSeq &seq) {
        constexpr std::size_t p = ((w - 1) / 32) + 1;

        std::array<std::uint_least32_t, n * p> a;
        seq.generate(a.begin(), a.end());

        auto iter = a.begin();
        for (std::size_t n_idx = 0; n_idx < n; n_idx++) {
            UIntType val = 0;
            for (std::size_t p_idx = 0; p_idx < p; ++p_idx) {
                val += static_cast<UIntType>(*iter++) << 32 * p_idx;
            }
            this->K[n_idx] = val & mask_result;
        }
        this->init();
    }

    void set_counter(const std::array<UIntType, n> &c) {
        for (std::size_t i = 0; i < n; ++i) {
            this->X[i] = c[i] & mask_result;
        }
        this->j = n - 1;
        this->operator++();
    }

    inline auto operator*() const -> value_type { return this->Y[this->j]; };

    inline auto operator++() -> threefry_engine & {
        ++this->j;
        if (this->j == n) {
            this->generate();
            this->increase_counter();
            this->j = 0;
        }
        return *this;
    }

    inline auto operator--() -> threefry_engine & {
        if (this->j == 0) {
            this->j = n;
            this->decrease_counter();
            this->decrease_counter();
            this->generate();
            this->increase_counter();
        }
        --this->j;
        return *this;
    }

    using common_engine<UIntType, w, threefry_engine>::operator++;
    using common_engine<UIntType, w, threefry_engine>::operator--;
    using common_engine<UIntType, w, threefry_engine>::operator==;

    auto operator==(const threefry_engine &rhs) const -> bool {
        return (this->X == rhs.X) && (this->K == rhs.K) && (this->Y == rhs.Y) &&
               (this->j == rhs.j);
    }

  private:
    std::array<UIntType, n> X{0};
    std::array<UIntType, n + 1> K{};
    std::array<UIntType, n> Y{0};
    std::size_t j{n - 1};
};

using threefry2x32 = threefry_engine<std::uint_fast32_t, 32, 2, 20, 0x1BD11BDA,
                                     13, 15, 26, 6, 17, 29, 16, 24>;

using threefry4x32 =
    threefry_engine<std::uint_fast32_t, 32, 4, 20, 0x1BD11BDA, 10, 26, 11, 21,
                    13, 27, 23, 5, 6, 20, 17, 11, 25, 10, 18, 20>;

using threefry2x64 =
    threefry_engine<std::uint_fast64_t, 64, 2, 20, 0x1BD11BDAA9FC1A22, 16, 42,
                    12, 31, 16, 32, 24, 21>;

using threefry4x64 =
    threefry_engine<std::uint_fast64_t, 64, 4, 20, 0x1BD11BDAA9FC1A22, 14, 16,
                    52, 57, 23, 40, 5, 37, 25, 33, 46, 12, 58, 22, 32, 32>;

} // namespace adhoc

#endif // BRNG_THREEFRY_ENGINE
