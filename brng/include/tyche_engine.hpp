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

#ifndef BRNG_TYCHE_ENGINE
#define BRNG_TYCHE_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/mask.hpp"
#include "tools/seed_seq_filler.hpp"

#include <cstddef>
#include <cstdint>

namespace adhoc {

template <class UIntType, std::size_t w, bool open_rand = false>
class tyche_engine final
    : public common_engine<UIntType, w, tyche_engine<UIntType, w, open_rand>> {
  public:
    static constexpr std::size_t word_size = w;
    using value_type = UIntType;
    static constexpr auto mask_result = tyche_engine::max();
    static constexpr UIntType default_seed = 5489U;

    tyche_engine() : tyche_engine(default_seed) {}

    explicit tyche_engine(std::uint64_t seed, std::uint32_t idx = 0U) {
        constexpr uint32_t global_seed = 0xAAAAAAAA;

        if constexpr (open_rand) {
            seed ^= global_seed;
        }
        constexpr auto mask32 = mask<std::uint64_t>(32);
        this->a = static_cast<uint32_t>(seed >> 32U);
        this->b = static_cast<uint32_t>(seed & mask32);
        this->c = 2654435769U;
        this->d = 1367130551U;
        this->d = this->d ^ idx;

        this->discard(20);
        this->operator++();
    }

    template <class SeedSeq> explicit tyche_engine(SeedSeq &seq) {
        seed_seq_filler<w, 4> seq_filler(seq);
        this->a = seq_filler.template get<UIntType>();
        this->b = seq_filler.template get<UIntType>();
        this->c = seq_filler.template get<UIntType>();
        this->d = seq_filler.template get<UIntType>();
        this->operator++();
    }

    inline auto operator*() const -> value_type {
        if constexpr (open_rand) {
            return this->a;
        } else {
            return this->b;
        }
    };

    inline auto operator++() -> tyche_engine & {
        if constexpr (open_rand) {
            this->b = rotl<word_size>(this->b, 7) ^ this->c;
            this->c -= this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->d = rotl<word_size>(this->d, 8) ^ this->a;
            this->a -= this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
            this->b = rotl<word_size>(this->b, 12) ^ this->c;
            this->c -= this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->d = rotl<word_size>(this->d, 16) ^ this->a;
            this->a -= this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
        } else {
            this->a += this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
            this->d = rotl<word_size>(this->d ^ this->a, 16);
            this->c += this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->b = rotl<word_size>(this->b ^ this->c, 12);
            this->a += this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
            this->d = rotl<word_size>(this->d ^ this->a, 8);
            this->c += this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->b = rotl<word_size>(this->b ^ this->c, 7);
        }
        return *this;
    }

    inline auto operator--() -> tyche_engine & {
        if constexpr (open_rand) {
            this->a += this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
            this->d = rotr<word_size>(this->d ^ this->a, 16);
            this->c += this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->b = rotr<word_size>(this->b ^ this->c, 12);
            this->a += this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
            this->d = rotr<word_size>(this->d ^ this->a, 8);
            this->c += this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->b = rotr<word_size>(this->b ^ this->c, 7);
        } else {
            this->b = rotr<word_size>(this->b, 7) ^ this->c;
            this->c -= this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->d = rotr<word_size>(this->d, 8) ^ this->a;
            this->a -= this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
            this->b = rotr<word_size>(this->b, 12) ^ this->c;
            this->c -= this->d;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            this->d = rotr<word_size>(this->d, 16) ^ this->a;
            this->a -= this->b;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
        }
        return *this;
    }

    using common_engine<UIntType, w, tyche_engine>::operator++;
    using common_engine<UIntType, w, tyche_engine>::operator--;
    using common_engine<UIntType, w, tyche_engine>::operator==;

    auto operator==(const tyche_engine &rhs) const -> bool {
        return (this->a == rhs.a) && (this->b == rhs.b) && (this->c == rhs.c) &&
               (this->d == rhs.d);
    }

  private:
    UIntType a{0};
    UIntType b{0};
    UIntType c{0};
    UIntType d{0};
};

using tyche = tyche_engine<std::uint_fast32_t, 32>;
using tyche_openrand = tyche_engine<std::uint_fast32_t, 32, true>;

} // namespace adhoc

#endif // BRNG_TYCHE_ENGINE
