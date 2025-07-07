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

#ifndef BRNG_EFIIX_ENGINE
#define BRNG_EFIIX_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/seed_seq_filler.hpp"

#include "jsf_engine.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, int rot_amount>
class efiix_engine final
    : public common_engine<UIntType, w, efiix_engine<UIntType, w, rot_amount>> {
  public:
    using value_type = UIntType;

    static constexpr UIntType default_seed = 5489U;
    static constexpr UIntType mask_result = efiix_engine::max();

    efiix_engine() : efiix_engine(default_seed) {}

    explicit efiix_engine(std::uint64_t s1, std::uint64_t s2 = 0U,
                          std::uint64_t s3 = 0U, std::uint64_t s4 = 0U) {
        s2 = s1;
        s3 = s1;
        s4 = s1;
        if constexpr (w == 8) {
            this->a = this->b = this->c = this->i = 0;
            for (int x = 0; x < 1 + (64 / 8); x++) {
                this->a += UIntType(s1);
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    this->a &= mask_result;
                }
                this->b += UIntType(s2);
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    this->b &= mask_result;
                }
                this->c += UIntType(s3);
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    this->c &= mask_result;
                }
                this->iteration_table[0] += UIntType(s4);
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    this->iteration_table[0] &= mask_result;
                }
                s1 >>= 8U;
                s2 >>= 8U;
                s3 >>= 1U;
                s4 >>= 1U;
                enum {
                    LESSER_SIZE = ITERATION_SIZE > INDIRECTION_SIZE
                                      ? INDIRECTION_SIZE
                                      : ITERATION_SIZE
                };
                enum {
                    GREATER_SIZE = ITERATION_SIZE > INDIRECTION_SIZE
                                       ? ITERATION_SIZE
                                       : INDIRECTION_SIZE
                };
                if (!x) {
                    unsigned long mask = 1;
                    while (mask < LESSER_SIZE - 1) {
                        for (unsigned long y = 0; y <= mask; y++) {
                            UIntType iterated =
                                this->iteration_table[this->i & mask &
                                                      (ITERATION_SIZE - 1)];
                            UIntType indirect =
                                this->indirection_table[this->c & mask &
                                                        (INDIRECTION_SIZE - 1)];
                            this->indirection_table[c & mask] = iterated + a;
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->indirection_table[c & mask] &=
                                    mask_result;
                            }
                            this->iteration_table[i & mask] = indirect;
                            UIntType const old = this->a ^ this->b;
                            this->a = this->b + this->i++;
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->a &= mask_result;
                            }
                            this->b = this->c + indirect;
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->b &= mask_result;
                            }
                            this->c = old + rotl<w>(this->c, 3);
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->c &= mask_result;
                            }
                            this->indirection_table[y + mask + 1] =
                                this->b ^ iterated;
                        }
                        mask = (mask << 1U) | 1U;
                    }
                    while (mask < GREATER_SIZE - 1) {
                        for (unsigned long y = 0; y <= mask; y++) {
                            UIntType iterated =
                                this->iteration_table[this->i & mask &
                                                      (ITERATION_SIZE - 1)];
                            UIntType indirect =
                                this->indirection_table[c & mask &
                                                        (INDIRECTION_SIZE - 1)];
                            this->indirection_table[this->c & mask] =
                                iterated + this->a;
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->indirection_table[this->c & mask] &=
                                    mask_result;
                            }
                            this->iteration_table[this->i & mask] = indirect;
                            UIntType const old = this->a ^ this->b;
                            this->a = this->b + this->i++;
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->a &= mask_result;
                            }
                            this->b = this->c + indirect;
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->b &= mask_result;
                            }
                            this->c = old + rotl<w>(this->c, 3);
                            if constexpr (w < std::numeric_limits<
                                                  UIntType>::digits) {
                                this->c &= mask_result;
                            }
                        }
                        mask = (mask << 1U) | 1U;
                    }
                } else {
                    unsigned long y = 0;
                    for (; y < 20; y++) {
                        this->operator++();
                    }
                    for (; y < GREATER_SIZE; y++) {
                        this->operator++();
                    }
                }
            }
        } else {
            arbee seeder(s1, s2, s3, s4, 12);

            for (auto &indirection_val : this->indirection_table) {
                indirection_val = static_cast<UIntType>(seeder());
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    indirection_val &= mask_result;
                }
            }

            this->i = static_cast<UIntType>(seeder());
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->i &= mask_result;
            }
            for (std::size_t j = 0; j < ITERATION_SIZE; j++) {
                this->iteration_table[(j + this->i) % ITERATION_SIZE] =
                    static_cast<UIntType>(seeder());
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    this->iteration_table[(j + this->i) % ITERATION_SIZE] &=
                        mask_result;
                }
            }

            this->a = static_cast<UIntType>(seeder());
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->a &= mask_result;
            }
            this->b = static_cast<UIntType>(seeder());
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->b &= mask_result;
            }
            this->c = static_cast<UIntType>(seeder());
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->c &= mask_result;
            }
            for (std::size_t j = 0; j < 64; j++) {
                this->operator++();
            }

            seeder();
            s1 += seeder();
            s2 += seeder();
            s3 += seeder();
            arbee seeder2(s1 ^ this->a, s2 ^ this->b, s3 ^ this->c, ~s4, 12);
            for (auto &indirection_val : this->indirection_table) {
                indirection_val ^= static_cast<UIntType>(seeder2());
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    indirection_val &= mask_result;
                }
            }

            for (std::size_t j = 0; j < ITERATION_SIZE + 16; j++) {
                this->operator++();
            }
        }
    }

    template <class SeedSeq> explicit efiix_engine(SeedSeq &seq) {
        seed_seq_filler<w, 4 + ITERATION_SIZE + INDIRECTION_SIZE> seq_filler(
            seq);

        this->a = seq_filler.template get<UIntType>();
        this->b = seq_filler.template get<UIntType>();
        this->c = seq_filler.template get<UIntType>();
        this->i = seq_filler.template get<UIntType>();
        seq_filler.template generate<UIntType>(this->iteration_table.begin(),
                                               this->iteration_table.end());
        seq_filler.template generate<UIntType>(this->indirection_table.begin(),
                                               this->indirection_table.end());
    }

    inline auto operator*() const -> value_type {
        UIntType const iterated =
            this->iteration_table[static_cast<std::size_t>(this->i) %
                                  ITERATION_SIZE];
        UIntType const indirect =
            this->indirection_table[static_cast<std::size_t>(this->c) %
                                    INDIRECTION_SIZE];
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            return ((this->c + indirect) ^ iterated) & mask_result;
        } else {
            return (this->c + indirect) ^ iterated;
        }
    }

    inline auto operator++() -> efiix_engine & {
        UIntType const iterated =
            this->iteration_table[static_cast<std::size_t>(this->i) %
                                  ITERATION_SIZE];
        UIntType const indirect =
            this->indirection_table[static_cast<std::size_t>(this->c) %
                                    INDIRECTION_SIZE];
        this->indirection_table[static_cast<std::size_t>(this->c) %
                                INDIRECTION_SIZE] = iterated + this->a;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->indirection_table[static_cast<std::size_t>(this->c) %
                                    INDIRECTION_SIZE] &= mask_result;
        }
        this->iteration_table[static_cast<std::size_t>(this->i) %
                              ITERATION_SIZE] = indirect;
        UIntType const old = this->a ^ this->b;
        this->a = this->b + this->i++;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a &= mask_result;
        }
        this->b = this->c + indirect;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b &= mask_result;
        }
        this->c = old + rotl<w>(this->c, rot_amount);
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->c &= mask_result;
        }
        return *this;
    }

    inline auto operator--() -> efiix_engine & {
        --this->i;
        UIntType const indirect =
            this->iteration_table[this->i % ITERATION_SIZE];
        UIntType const new_c = this->c;
        UIntType const new_b = this->b;
        this->c = this->b - indirect;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->c &= mask_result;
        }
        UIntType const old = new_c - rotl<w>(this->c, rot_amount);
        this->b = this->a - this->i;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b &= mask_result;
        }
        this->a = old ^ this->b;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a &= mask_result;
        }
        UIntType iterated =
            this->indirection_table[this->c % INDIRECTION_SIZE] - this->a;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            iterated &= mask_result;
        }
        this->indirection_table[this->c % INDIRECTION_SIZE] = indirect;
        this->iteration_table[this->i % ITERATION_SIZE] = iterated;
        return *this;
    }

    using common_engine<UIntType, w, efiix_engine>::operator++;
    using common_engine<UIntType, w, efiix_engine>::operator--;
    using common_engine<UIntType, w, efiix_engine>::operator==;

    auto operator==(const efiix_engine &rhs) const -> bool {
        return (this->iteration_table == rhs.iteration_table) &&
               (this->indirection_table == rhs.indirection_table) &&
               (this->i == rhs.i) && (this->a == rhs.a) && (this->b == rhs.b) &&
               (this->c == rhs.c);
    }

  private:
    static constexpr std::size_t ITERATION_SIZE_L2 = 5;
    static constexpr std::size_t ITERATION_SIZE = 1U << ITERATION_SIZE_L2;
    static constexpr std::size_t INDIRECTION_SIZE_L2 = 4;
    static constexpr std::size_t INDIRECTION_SIZE = 1U << INDIRECTION_SIZE_L2;

    std::array<UIntType, ITERATION_SIZE> iteration_table{};
    std::array<UIntType, INDIRECTION_SIZE> indirection_table{};
    UIntType i, a, b, c;
};

using efiix64x48 = efiix_engine<std::uint64_t, 64, 25>;
using efiix32x48 = efiix_engine<std::uint32_t, 32, 13>;
using efiix16x48 = efiix_engine<std::uint16_t, 16, 7>;
using efiix8x48 = efiix_engine<std::uint8_t, 8, 3>;

} // namespace adhoc

#endif // BRNG_EFIIX_ENGINE
