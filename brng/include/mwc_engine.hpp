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

#ifndef BRNG_MWC
#define BRNG_MWC

#include "tools/circular_buffer.hpp"
#include "tools/common_engine.hpp"
#include "tools/uint128.hpp"

#include <cstdint>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t bitlength, UIntType m>
class mcw_engine final
    : public common_engine<UIntType, w, mcw_engine<UIntType, w, bitlength, m>> {
  public:
    static_assert(bitlength % w == 0);
    static constexpr std::size_t size = bitlength / w;
    static constexpr UIntType mask_result = mcw_engine::max();

    using upgraded_type = std::conditional_t<
        w <= 8, std::uint_fast16_t,
        std::conditional_t<
            w <= 16, std::uint_fast32_t,
            std::conditional_t<w <= 32, std::uint_fast64_t, uint128>>>;
    using value_type = UIntType;

    template <class... value_type> explicit mcw_engine(value_type... values) {
        static_assert(sizeof...(values) <= size);

        auto make_array =
            []<class... Args>(Args... args) -> std::array<UIntType, size> {
            return {static_cast<UIntType>(args)...};
        };
        auto &data = state.data();
        data = make_array(values...);
        data.back() %= m;
    }

    inline auto operator*() const -> value_type { return this->state.at(); }

    inline auto operator++() -> mcw_engine & {
        --this->state;
        const upgraded_type temp =
            (static_cast<upgraded_type>(this->state.template at<size - 1>()) *
             m) +
            this->state.at();

        if constexpr (w < std::numeric_limits<UIntType>::digits) {
            this->state.at() = temp & mask_result;
        } else {
            this->state.at() = temp;
        }
        this->state.template at<size - 1>() = temp >> w;

        return *this;
    }

    inline auto operator--() -> mcw_engine & {
        // this is a little tricky so worth a comment:
        // MCW is equivalent to LCG with the same multiplier and modulo
        // (a*(2^64) - 1) so naturally the inverse would be using 2^64 as
        // multiplier (a*(2^64) = 1[mod(a*(2^64) - 1)]) and using the MCW<->
        // LCG equivalence again on (2^64) as multiplier and same modulo
        // (a*(2^64) - 1), we have this MCW code:
        constexpr upgraded_type m_inv = static_cast<upgraded_type>(1) << w;
        const upgraded_type temp =
            (static_cast<upgraded_type>(this->state.template at<size - 1>()) *
             m_inv) +
            this->state.at();

        this->state.template at<size - 1>() = temp / m;
        this->state.at() = temp % m;
        ++this->state;

        return *this;
    }

    using common_engine<UIntType, w, mcw_engine>::operator++;
    using common_engine<UIntType, w, mcw_engine>::operator--;
    using common_engine<UIntType, w, mcw_engine>::operator==;

    auto operator==(const mcw_engine &rhs) const -> bool {
        return this->state == rhs.state;
    }

  private:
    circular_buffer<UIntType, size> state{};
};

using mcw128 = mcw_engine<std::uint_fast64_t, 64, 128, 0xffebb71d94fcdaf9>;
using mcw192 = mcw_engine<std::uint_fast64_t, 64, 192, 0xffa04e67b3c95d86>;
using mcw256 = mcw_engine<std::uint_fast64_t, 64, 256, 0xfff62cf2ccc0cdaf>;

} // namespace adhoc

#endif // BRNG_MRG31k3p
