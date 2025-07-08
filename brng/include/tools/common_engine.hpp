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

#ifndef BRNG_TOOLS_COMMON_ENGINE
#define BRNG_TOOLS_COMMON_ENGINE

#include "tools/mask.hpp"

#include <cstddef>
#include <iterator>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, class RNG> class common_engine {
  public:
    static constexpr std::size_t word_size = w;
    static_assert(word_size <= std::numeric_limits<UIntType>::digits);
    static_assert(std::is_unsigned_v<UIntType>);

    using difference_type = std::ptrdiff_t;
    using value_type = UIntType;
    using result_type = UIntType;

    inline auto operator*() const -> value_type {
        return (static_cast<RNG const *>(this))->operator*();
    }

    inline auto operator++() -> RNG & {
        return (static_cast<RNG *>(this))->operator++();
    }

    inline auto operator--() -> RNG & {
        return (static_cast<RNG *>(this))->operator--();
    }

    // creates a copy. should be avoided!
    inline auto operator++(int) -> RNG {
        auto const tmp = (*static_cast<RNG *>(this));
        ++*this;
        return tmp;
    }

    // creates a copy. should be avoided!
    inline auto operator--(int) -> RNG {
        auto const tmp = (*static_cast<RNG *>(this));
        --*this;
        return tmp;
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            auto const result = this->operator*();
            this->operator++();
            return result;
        } else {
            return *(this->operator--());
        }
    }

    void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator++();
        }
    }

    static constexpr auto min() -> UIntType {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> UIntType {
        return mask<UIntType>(word_size);
    }

    auto operator!=(const common_engine &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

    constexpr auto operator==(std::default_sentinel_t /*unused*/) const
        -> bool {
        return false;
    }
};

} // namespace adhoc

#endif // BRNG_TOOLS_COMMON_ENGINE
