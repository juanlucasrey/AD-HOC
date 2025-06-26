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

#ifndef BRNG_JSF_ENGINE
#define BRNG_JSF_ENGINE

#include "tools/bit.hpp"
#include "tools/mask.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t p, std::size_t q,
          std::size_t r, bool arbee = false>
class jsf_engine final {
  public:
    using result_type = UIntType;
    static constexpr std::size_t word_size = w;

    static constexpr auto min() -> result_type {
        return static_cast<result_type>(0U);
    };

    static constexpr auto max() -> result_type {
        return mask<UIntType>(word_size);
    };

    static_assert(w <= std::numeric_limits<UIntType>::digits);
    static_assert(p < w);
    static_assert(q < w);
    static_assert(r < w);

    static constexpr UIntType m = jsf_engine::max();
    static constexpr UIntType default_seed_1 = 0xf1ea5eed;
    static constexpr UIntType default_seed_2 = 0xcafe5eed00000001ULL;
    static constexpr std::size_t default_seed_3 = 20U;

    jsf_engine()
        : jsf_engine(default_seed_1, default_seed_2, default_seed_2,
                     default_seed_2, default_seed_3) {}

    explicit jsf_engine(result_type seed1, result_type seed2 = default_seed_2,
                        result_type seed3 = default_seed_2,
                        result_type seed4 = default_seed_2,
                        std::size_t iters = default_seed_3) {
        this->state[0] = seed1 & m;
        this->state[1] = seed2 & m;
        this->state[2] = seed3 & m;
        this->state[3] = seed4 & m;
        if constexpr (arbee) {
            this->state[4] = 1;
        }
        for (std::size_t j = 0; j < iters; ++j) {
            this->operator()<true>();
        }
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            UIntType e = this->state[0];
            if constexpr (arbee) {
                e += rotl<w>(this->state[1], p);
            } else {
                e -= rotl<w>(this->state[1], p);
            }

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                e &= m;
            }

            this->state[0] = this->state[1] ^ rotl<w>(this->state[2], q);
            if constexpr (r) {
                this->state[1] = this->state[2] + rotl<w>(this->state[3], r);
            } else {
                this->state[1] = this->state[2] + this->state[3];
            }

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->state[1] &= m;
            }

            if constexpr (arbee) {
                this->state[2] = this->state[3] + e + this->state[4]++;
            } else {
                this->state[2] = this->state[3] + e;
            }

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->state[2] &= m;
            }
            this->state[3] = e + this->state[0];

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->state[3] &= m;
            }

            return this->state[3];
        } else {
            auto const result = this->state[3];
            UIntType e = this->state[3] - this->state[0];

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                e &= m;
            }

            if constexpr (arbee) {
                this->state[3] = this->state[2] - e - --this->state[4];
            } else {
                this->state[3] = this->state[2] - e;
            }

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->state[3] &= m;
            }

            if constexpr (r) {
                this->state[2] = this->state[1] - rotl<w>(this->state[3], r);
            } else {
                this->state[2] = this->state[1] - this->state[3];
            }

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->state[2] &= m;
            }

            this->state[1] = this->state[0] ^ rotl<w>(this->state[2], q);

            if constexpr (arbee) {
                this->state[0] = e - rotl<w>(this->state[1], p);
            } else {
                this->state[0] = e + rotl<w>(this->state[1], p);
            }

            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->state[0] &= m;
            }

            return result;
        }
    }

    template <bool FwdDirection = true> void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }

    auto operator==(const jsf_engine &rhs) const -> bool {
        return (this->state == rhs.state);
    }

    auto operator!=(const jsf_engine &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    std::array<UIntType, arbee ? 5 : 4> state{};
};

using jsf32na = jsf_engine<std::uint32_t, 32, 27, 17, 0>;
using jsf32nb = jsf_engine<std::uint32_t, 32, 9, 16, 0>;
using jsf32nc = jsf_engine<std::uint32_t, 32, 9, 24, 0>;
using jsf32nd = jsf_engine<std::uint32_t, 32, 10, 16, 0>;
using jsf32ne = jsf_engine<std::uint32_t, 32, 10, 24, 0>;
using jsf32nf = jsf_engine<std::uint32_t, 32, 11, 16, 0>;
using jsf32ng = jsf_engine<std::uint32_t, 32, 11, 24, 0>;
using jsf32nh = jsf_engine<std::uint32_t, 32, 25, 8, 0>;
using jsf32ni = jsf_engine<std::uint32_t, 32, 25, 16, 0>;
using jsf32nj = jsf_engine<std::uint32_t, 32, 26, 8, 0>;
using jsf32nk = jsf_engine<std::uint32_t, 32, 26, 16, 0>;
using jsf32nl = jsf_engine<std::uint32_t, 32, 26, 17, 0>;
using jsf32nm = jsf_engine<std::uint32_t, 32, 27, 16, 0>;
using jsf32ra = jsf_engine<std::uint32_t, 32, 3, 14, 24>;
using jsf32rb = jsf_engine<std::uint32_t, 32, 3, 25, 15>;
using jsf32rc = jsf_engine<std::uint32_t, 32, 4, 15, 24>;
using jsf32rd = jsf_engine<std::uint32_t, 32, 6, 16, 28>;
using jsf32re = jsf_engine<std::uint32_t, 32, 7, 16, 27>;
using jsf32rf = jsf_engine<std::uint32_t, 32, 8, 14, 3>;
using jsf32rg = jsf_engine<std::uint32_t, 32, 11, 16, 23>;
using jsf32rh = jsf_engine<std::uint32_t, 32, 12, 16, 22>;
using jsf32ri = jsf_engine<std::uint32_t, 32, 12, 17, 23>;
using jsf32rj = jsf_engine<std::uint32_t, 32, 13, 16, 22>;
using jsf32rk = jsf_engine<std::uint32_t, 32, 15, 25, 3>;
using jsf32rl = jsf_engine<std::uint32_t, 32, 16, 9, 3>;
using jsf32rm = jsf_engine<std::uint32_t, 32, 17, 9, 3>;
using jsf32rn = jsf_engine<std::uint32_t, 32, 17, 27, 7>;
using jsf32ro = jsf_engine<std::uint32_t, 32, 19, 7, 3>;
using jsf32rp = jsf_engine<std::uint32_t, 32, 23, 15, 11>;
using jsf32rq = jsf_engine<std::uint32_t, 32, 23, 16, 11>;
using jsf32rr = jsf_engine<std::uint32_t, 32, 23, 17, 11>;
using jsf32rs = jsf_engine<std::uint32_t, 32, 24, 3, 16>;
using jsf32rt = jsf_engine<std::uint32_t, 32, 24, 4, 16>;
using jsf32ru = jsf_engine<std::uint32_t, 32, 25, 14, 3>;
using jsf32rv = jsf_engine<std::uint32_t, 32, 27, 16, 6>;
using jsf32rw = jsf_engine<std::uint32_t, 32, 27, 16, 7>;
using jsf32n = jsf32na;
using jsf32r = jsf32rq;
using jsf32 = jsf32n;
using jsf64na = jsf_engine<std::uint64_t, 64, 39, 11, 0>;
using jsf64ra = jsf_engine<std::uint64_t, 64, 7, 13, 37>;
using jsf64n = jsf64na;
using jsf64r = jsf64ra;
using jsf64 = jsf64r;
using jsf16na = jsf_engine<std::uint64_t, 64, 13, 8, 0>;
using jsf16 = jsf16na;
using jsf8na = jsf_engine<std::uint8_t, 8, 1, 4, 0>;
using jsf8 = jsf8na;
using arbee = adhoc::jsf_engine<std::uint64_t, 64, 45, 13, 37, true>;

} // namespace adhoc

#endif // BRNG_JSF_ENGINE
