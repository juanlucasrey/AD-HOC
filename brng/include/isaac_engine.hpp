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

#ifndef BRNG_ISAAC_ENGINE
#define BRNG_ISAAC_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/mask.hpp"
#include "tools/seed_seq_filler.hpp"
#include "tools/unshift.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, int shift1, int shift2, int shift3,
          int shift4, int shift5, bool flip, UIntType ratio>
class isaac_engine final
    : public common_engine<UIntType, w,
                           isaac_engine<UIntType, w, shift1, shift2, shift3,
                                        shift4, shift5, flip, ratio>> {
  private:
    static_assert(w == 64 || w == 32);
    static constexpr std::size_t power2_size = 8;
    static constexpr std::size_t size = 1U << power2_size;
    static constexpr std::size_t mask_size = mask<std::uint32_t>(power2_size);

  public:
    using value_type = UIntType;
    static constexpr auto mask_result = isaac_engine::max();
    static constexpr UIntType default_seed = 5489U;

    isaac_engine() : isaac_engine(default_seed) {}

    explicit isaac_engine(UIntType s1, UIntType s2 = 0U) {
        this->state[0] = s1;
        this->state[1] = s2;

        std::array<UIntType, 8> t{};
        t.fill(ratio);

        auto mix = [&t]() {
            if constexpr (w == 32) {
                t[0] ^= t[1] << 11;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[0] &= mask_result;
                }
                t[3] += t[0];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[3] &= mask_result;
                }
                t[1] += t[2];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[1] &= mask_result;
                }
                t[1] ^= t[2] >> 2;
                t[4] += t[1];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[4] &= mask_result;
                }
                t[2] += t[3];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[2] &= mask_result;
                }
                t[2] ^= t[3] << 8;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[2] &= mask_result;
                }
                t[5] += t[2];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[5] &= mask_result;
                }
                t[3] += t[4];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[3] &= mask_result;
                }
                t[3] ^= t[4] >> 16;
                t[6] += t[3];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[6] &= mask_result;
                }
                t[4] += t[5];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[4] &= mask_result;
                }
                t[4] ^= t[5] << 10;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[4] &= mask_result;
                }
                t[7] += t[4];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[7] &= mask_result;
                }
                t[5] += t[6];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[5] &= mask_result;
                }
                t[5] ^= t[6] >> 4;
                t[0] += t[5];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[0] &= mask_result;
                }
                t[6] += t[7];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[6] &= mask_result;
                }
                t[6] ^= t[7] << 8;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[6] &= mask_result;
                }
                t[1] += t[6];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[1] &= mask_result;
                }
                t[7] += t[0];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[7] &= mask_result;
                }
                t[7] ^= t[0] >> 9;
                t[2] += t[7];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[2] &= mask_result;
                }
                t[0] += t[1];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[0] &= mask_result;
                }
            } else {
                t[0] -= t[4];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[0] &= mask_result;
                }
                t[5] ^= t[7] >> 9;
                t[7] += t[0];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[7] &= mask_result;
                }
                t[1] -= t[5];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[1] &= mask_result;
                }
                t[6] ^= t[0] << 9;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[6] &= mask_result;
                }
                t[0] += t[1];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[0] &= mask_result;
                }
                t[2] -= t[6];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[2] &= mask_result;
                }
                t[7] ^= t[1] >> 23;
                t[1] += t[2];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[1] &= mask_result;
                }
                t[3] -= t[7];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[3] &= mask_result;
                }
                t[0] ^= t[2] << 15;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[0] &= mask_result;
                }
                t[2] += t[3];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[2] &= mask_result;
                }
                t[4] -= t[0];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[4] &= mask_result;
                }
                t[1] ^= t[3] >> 14;
                t[3] += t[4];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[3] &= mask_result;
                }
                t[5] -= t[1];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[5] &= mask_result;
                }
                t[2] ^= t[4] << 20;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[2] &= mask_result;
                }
                t[4] += t[5];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[4] &= mask_result;
                }
                t[6] -= t[2];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[6] &= mask_result;
                }
                t[3] ^= t[5] >> 17;
                t[5] += t[6];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[5] &= mask_result;
                }
                t[7] -= t[3];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[7] &= mask_result;
                }
                t[4] ^= t[6] << 14;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[4] &= mask_result;
                }
                t[6] += t[7];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    t[6] &= mask_result;
                }
            }
        };

        for (std::size_t i = 0; i < 4; ++i) {
            mix();
        }

        for (std::size_t pass = 0; pass < 2; pass++) {
            for (std::size_t i = 0; i < size; i += 8) {
                for (std::size_t x = 0; x < 8; x++) {
                    t[x] += this->state[i + x];
                    if constexpr (w < std::numeric_limits<UIntType>::digits) {
                        t[x] &= mask_result;
                    }
                }

                mix();

                for (std::size_t x = 0; x < 8; x++) {
                    this->state[i + x] = t[x];
                }
            }
        }
        this->operator++();
    }

    template <class SeedSeq> explicit isaac_engine(SeedSeq &seq) {
        seed_seq_filler<w, size + 3> seq_filler(seq);
        seq_filler.template generate<UIntType>(state.begin(), state.end());
        this->a = seq_filler.template get<UIntType>();
        this->b = seq_filler.template get<UIntType>();
        this->c = seq_filler.template get<UIntType>();
        this->operator++();
    }

    inline auto operator*() const -> value_type { return this->b; }

    inline auto operator++() -> isaac_engine & {
        constexpr std::array<int, 4> shift_vals = {shift1, shift2, shift3,
                                                   shift4};

        ++this->idx;
        this->idx &= mask_size;

        if (this->idx == 0) {
            this->b += ++this->c;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->b &= mask_result;
                this->c &= mask_result;
            }
        }

        UIntType const x = this->state[this->idx];

        this->a = (this->a ^ shift(this->a, shift_vals[this->idx % 4]));
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a &= mask_result;
        }
        if constexpr (flip) {
            if ((this->idx % 4) == 0) {
                this->a = ~this->a;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    this->a &= mask_result;
                }
            }
        }

        this->a += this->state[(this->idx + size / 2) & mask_size];
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a &= mask_result;
        }

        this->state[this->idx] =
            this->state[(x >> shift5) & mask_size] + this->a + this->b;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->state[this->idx] &= mask_result;
        }
        this->b =
            this->state[(this->state[this->idx] >> (8 + shift5)) & mask_size] +
            x;
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b &= mask_result;
        }

        return *this;
    }

    inline auto operator--() -> isaac_engine & {
        constexpr std::array<int, 4> shift_vals = {shift1, shift2, shift3,
                                                   shift4};
        UIntType x =
            this->b -
            this->state[(this->state[this->idx] >> (8 + shift5)) & mask_size];
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            x &= mask_result;
        }

        if (static_cast<std::size_t>((x >> shift5) & mask_size) == this->idx) {
            this->b = this->state[this->idx] - x - this->a;
        } else {
            this->b = this->state[this->idx] -
                      this->state[(x >> shift5) & mask_size] - this->a;
        }

        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->b &= mask_result;
        }

        this->a -= this->state[(this->idx + size / 2) & mask_size];
        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->a &= mask_result;
        }

        auto const shift_idx = this->idx % 4;

        if (shift_idx == 0) {
            if constexpr (flip) {
                this->a = ~this->a;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    this->a &= mask_result;
                }
            }
            this->a = unshift_xor<w, shift_vals[0]>(this->a);
        } else if (shift_idx == 1) {
            this->a = unshift_xor<w, shift_vals[1]>(this->a);
        } else if (shift_idx == 2) {
            this->a = unshift_xor<w, shift_vals[2]>(this->a);
        } else if (shift_idx == 3) {
            this->a = unshift_xor<w, shift_vals[3]>(this->a);
        }

        this->state[this->idx] = x;

        if (this->idx == 0) {
            this->b -= this->c;
            --this->c;
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->b &= mask_result;
                this->c &= mask_result;
            }
        }

        --this->idx;
        this->idx &= mask_size;

        return *this;
    }

    using common_engine<UIntType, w, isaac_engine>::operator++;
    using common_engine<UIntType, w, isaac_engine>::operator--;
    using common_engine<UIntType, w, isaac_engine>::operator==;

    auto operator==(const isaac_engine &rhs) const -> bool {
        return (this->a == rhs.a) && (this->b == rhs.b) && (this->c == rhs.c) &&
               (this->idx == rhs.idx) && (this->state == rhs.state);
    }

  private:
    std::array<UIntType, size> state{};
    UIntType a{0};
    UIntType b{0};
    UIntType c{0};
    std::size_t idx = mask_size;
};

using isaac32x256 = isaac_engine<std::uint_fast32_t, 32, 13, -6, 2, -16, 2,
                                 false, 0x9e3779b9UL>;

using isaac64x256 = isaac_engine<std::uint_fast64_t, 64, 21, -5, 12, -33, 3,
                                 true, 0x9e3779b97f4a7c13ULL>;

} // namespace adhoc

#endif // BRNG_ISAAC_ENGINE
