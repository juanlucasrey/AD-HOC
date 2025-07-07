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

#ifndef BRNG_LXM_ENGINE
#define BRNG_LXM_ENGINE

#include "linear_congruential_engine.hpp"
#include "tools/common_engine.hpp"
#include "xoshiro_engine.hpp"

#include <cstddef>

namespace adhoc {

template <class UIntType, std::size_t w>
class lxm_engine final
    : public common_engine<UIntType, w, lxm_engine<UIntType, w>> {
  public:
    using value_type = UIntType;
    static constexpr UIntType mask_result = lxm_engine::max();
    static constexpr UIntType default_seed = 0x3EE517E67D1C29DEULL;

    lxm_engine() : lxm_engine(default_seed) {}

    explicit lxm_engine(UIntType lcg_seed,
                        UIntType xseed1 = 0x5C572D54940C542EULL,
                        UIntType xseed2 = 0x2710CDA3D55747B4ULL,
                        UIntType xseed3 = 0x631BD5B8C716B444ULL,
                        UIntType xseed4 = 0xD76F337C89CCCBAFULL)
        : lcg(lcg_seed), x(xseed1, xseed2, xseed3, xseed4) {
        --lcg;
    }

    template <class SeedSeq> explicit lxm_engine(SeedSeq &seq) {}

    inline auto operator*() const -> value_type {
        std::uint64_t key = x.get0() + *lcg;
        // splitmix64::operator*
        key = (key ^ (key >> 30U)) * 0xbf58476d1ce4e5b9;
        key = (key ^ (key >> 27U)) * 0x94d049bb133111eb;
        return key ^ (key >> 31U);
    }

    inline auto operator++() -> lxm_engine & {
        ++x;
        ++lcg;
        return *this;
    }

    inline auto operator--() -> lxm_engine & {
        --x;
        --lcg;
        return *this;
    }

    using common_engine<UIntType, w, lxm_engine>::operator++;
    using common_engine<UIntType, w, lxm_engine>::operator--;
    using common_engine<UIntType, w, lxm_engine>::operator==;

    auto operator==(const lxm_engine &rhs) const -> bool {
        return (this->lcg == rhs.lcg) && (this->x == rhs.x);
    }

  private:
    linear_congruential_engine<std::uint_fast64_t, 64, 2862933555777941757ULL,
                               3037000493ULL, 0>
        lcg;
    xoshiro256plus x;
};

} // namespace adhoc

#endif // BRNG_LXM_ENGINE
