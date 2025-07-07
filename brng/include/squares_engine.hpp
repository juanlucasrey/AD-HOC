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

#ifndef BRNG_SQUARES_ENGINE
#define BRNG_SQUARES_ENGINE

#include "tools/common_engine.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>

namespace adhoc {

template <class UIntType, std::size_t w>
class squares_engine final
    : public common_engine<UIntType, w, squares_engine<UIntType, w>> {

  public:
    using value_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::uint64_t default_seed_1 = 0xb504f333f9de6484UL;
    static constexpr UIntType mask_result = squares_engine::max();

    squares_engine() : squares_engine(default_seed_1) {}

    squares_engine(std::uint64_t key) : state(key) {}

    inline auto operator*() const -> value_type {
        std::uint64_t const y = ctr * this->state;
        std::uint64_t const z = y + this->state;
        std::uint64_t x = (y * y) + y;

        x = std::rotl(x, 32);
        x = x * x + z;
        x = std::rotl(x, 32);
        x = x * x + y;
        x = std::rotl(x, 32);
        if constexpr (w == 64) {
            x = x * x + z;
            std::uint64_t const t = x;
            x = std::rotl(x, 32);
            return t ^ ((x * x + y) >> 32U);
        } else {
            return (x * x + z) >> 32U;
        }
    }

    inline auto operator++() -> squares_engine & {
        ++this->ctr;
        return *this;
    }

    inline auto operator--() -> squares_engine & {
        --this->ctr;
        return *this;
    }

    using common_engine<UIntType, w, squares_engine>::operator++;
    using common_engine<UIntType, w, squares_engine>::operator--;
    using common_engine<UIntType, w, squares_engine>::operator==;

    auto operator==(const squares_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->ctr == rhs.ctr);
    }

  private:
    std::uint64_t state{};
    std::uint64_t ctr{};
};

using squares64 = squares_engine<std::uint_fast64_t, 64>;
using squares32 = squares_engine<std::uint_fast32_t, 32>;

} // namespace adhoc

#endif // BRNG_SQUARES_ENGINE
