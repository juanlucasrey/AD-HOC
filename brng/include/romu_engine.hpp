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

#ifndef BRNG_ROMU_ENGINE
#define BRNG_ROMU_ENGINE

#include "tools/common_engine.hpp"
#include "tools/modular_arithmetic.hpp"

#include <array>
#include <bit>
#include <cstddef>

namespace adhoc {

template <class UIntType, std::size_t w, class UIntType_int, std::size_t w_int,
          std::size_t s, UIntType_int mult, int rot1, int rot2 = 0>
class romu_engine final
    : public common_engine<
          UIntType, w,
          romu_engine<UIntType, w, UIntType_int, w_int, s, mult, rot1, rot2>> {
    static_assert(s <= 4);
    // only mono state is allowed to have different internal state
    static_assert(s == 1 || w == w_int);
    static_assert(rot1 < w_int);
    static_assert(rot2 < w_int);

  public:
    using value_type = UIntType;

    static constexpr UIntType_int default_seed = 0x9f57c403d06c42fcUL;

    romu_engine() : romu_engine(default_seed) {}

    explicit romu_engine(UIntType_int seed1, UIntType_int seed2 = 0,
                         UIntType_int seed3 = 0, UIntType_int seed4 = 0) {
        if constexpr (s == 1) {
            this->state[0] = (seed1 & 0x1fffffffU) + 1156979152U;
        } else if constexpr (s == 2) {
            this->state[0] = seed1;
            this->state[1] = seed2;
        } else if constexpr (s == 3) {
            this->state[0] = seed1;
            this->state[1] = seed2;
            this->state[2] = seed3;
        } else if constexpr (s == 4) {
            this->state[0] = seed1;
            this->state[1] = seed2;
            this->state[2] = seed3;
            this->state[3] = seed4;
        }
    }

    inline auto operator*() const -> value_type {
        if constexpr (s == 1) {
            return this->state[0] >> 16;
        } else if constexpr (s == 2 || s == 3) {
            return this->state[0];
        } else if constexpr (s == 4) {
            return this->state[1];
        }
    }

    inline auto operator++() -> romu_engine & {
        if constexpr (s == 1) {
            this->state[0] *= mult;
            this->state[0] = std::rotl(this->state[0], rot1);
        } else if constexpr (s == 2 && rot2) {
            UIntType const xp = this->state[0];
            this->state[0] = mult * this->state[1];
            this->state[1] = std::rotl(this->state[1], rot1) +
                             std::rotl(this->state[1], rot2) - xp;
        } else if constexpr (s == 2 && rot2 == 0) {
            UIntType const xp = this->state[0];
            this->state[0] = mult * this->state[1];
            this->state[1] = this->state[1] - xp;
            this->state[1] = std::rotl(this->state[1], rot1);
        } else if constexpr (s == 3) {
            UIntType const xp = this->state[0];
            UIntType const yp = this->state[1];
            UIntType const zp = this->state[2];
            this->state[0] = mult * zp;
            this->state[1] = yp - xp;
            this->state[1] = std::rotl(this->state[1], rot1);
            this->state[2] = zp - yp;
            this->state[2] = std::rotl(this->state[2], rot2);
        } else if constexpr (s == 4) {
            UIntType const wp = this->state[0];
            UIntType const xp = this->state[1];
            UIntType const yp = this->state[2];
            UIntType const zp = this->state[3];
            this->state[0] = mult * zp;
            this->state[1] = zp + std::rotl(wp, rot1);
            this->state[2] = yp - xp;
            this->state[3] = yp + wp;
            this->state[3] = std::rotl(this->state[3], rot2);
        }

        return *this;
    }

    inline auto operator--() -> romu_engine & {
        constexpr auto mult_inv =
            modular_multiplicative_inverse_pow2(w_int, mult);

        if constexpr (s == 1) {
            this->state[0] = std::rotr(this->state[0], rot1);
            this->state[0] *= mult_inv;
        } else if constexpr (s == 2 && rot2) {
            UIntType const yp = this->state[1];
            this->state[1] = mult_inv * this->state[0];
            this->state[0] = std::rotl(this->state[1], rot1) +
                             std::rotl(this->state[1], rot2) - yp;
        } else if constexpr (s == 2 && rot2 == 0) {
            UIntType const yp_mn_xp = std::rotr(this->state[1], rot1);
            this->state[1] = mult_inv * this->state[0];
            this->state[0] = this->state[1] - yp_mn_xp;
        } else if constexpr (s == 3) {
            UIntType const zp_mn_yp = std::rotr(this->state[2], rot2);
            this->state[2] = this->state[0] * mult_inv;
            UIntType const yp_mn_xp = std::rotr(this->state[1], rot1);
            this->state[1] = this->state[2] - zp_mn_yp;
            this->state[0] = this->state[1] - yp_mn_xp;
        } else if constexpr (s == 4) {
            UIntType const yp_pl_wp = std::rotr(this->state[3], rot2);
            this->state[3] = this->state[0] * mult_inv;
            this->state[0] = std::rotr(this->state[1] - this->state[3], rot1);
            UIntType const yp_mn_xp = this->state[2];
            this->state[2] = yp_pl_wp - this->state[0];
            this->state[1] = this->state[2] - yp_mn_xp;
        }
        return *this;
    }

    using common_engine<UIntType, w, romu_engine>::operator++;
    using common_engine<UIntType, w, romu_engine>::operator--;
    using common_engine<UIntType, w, romu_engine>::operator==;

    auto operator==(const romu_engine &rhs) const -> bool {
        return (this->state == rhs.state);
    }

  private:
    std::array<UIntType_int, s> state{};
};

using romuquad = romu_engine<std::uint64_t, 64, std::uint64_t, 64, 4,
                             15241094284759029579ULL, 52, 19>;

using romutrio = romu_engine<std::uint64_t, 64, std::uint64_t, 64, 3,
                             15241094284759029579ULL, 12, 44>;

using romuduo = romu_engine<std::uint64_t, 64, std::uint64_t, 64, 2,
                            15241094284759029579ULL, 36, 15>;

using romuduojr = romu_engine<std::uint64_t, 64, std::uint64_t, 64, 2,
                              15241094284759029579ULL, 27>;

using romuquad32 =
    romu_engine<std::uint32_t, 32, std::uint32_t, 32, 4, 3323815723U, 26, 9>;

using romutrio32 =
    romu_engine<std::uint32_t, 32, std::uint32_t, 32, 3, 3323815723U, 6, 22>;

using romumono32 =
    romu_engine<std::uint16_t, 16, std::uint32_t, 32, 1, 3611795771U, 12>;

} // namespace adhoc

#endif // BRNG_ROMU_ENGINE
