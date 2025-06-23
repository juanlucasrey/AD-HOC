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

#ifndef BRNG_TOOLS_SEED_SEQ_FILLER
#define BRNG_TOOLS_SEED_SEQ_FILLER

#include "mask.hpp"

#include <array>
#include <bit>
#include <cstdint>

namespace adhoc {

// this class adapts a SeedSeq outputting 32 bit words to any power of 2 sized
// words (8, 16, 32, 64, 128, ...)
template <std::size_t w, std::size_t n> class seed_seq_filler final {
    static_assert(std::has_single_bit(w));

  public:
    template <class SeedSeq> explicit seed_seq_filler(SeedSeq &seq) {
        seq.generate(seeds.begin(), seeds.end());
    }

    template <class UIntType> auto get() -> UIntType {
        static_assert(w <= std::numeric_limits<UIntType>::digits);

        UIntType result = 0;
        if constexpr (w == 32) {
            result = seeds.at(idx++);
        } else if constexpr (w < 32) {
            result = (seeds.at(idx) >> bits_used) & mask<UIntType>(w);
            bits_used += w;
            if (bits_used == 32) {
                bits_used = 0;
                ++idx;
            }
        } else {
            constexpr std::size_t ratio = w / 32;
            for (std::size_t i = 0; i < ratio; ++i) {
                result |= static_cast<UIntType>(seeds.at(idx++)) << (32 * i);
            }
        }
        return result;
    }

    template <class UIntType, typename RandomIt>
    void generate(RandomIt first, RandomIt last) {
        while (first != last) {
            *first = this->get<UIntType>();
            ++first;
        }
    }

  private:
    // add (32 - 1) to make sure we overshot
    std::array<std::uint32_t, (n * w + (32 - 1)) / 32> seeds{};
    std::size_t idx{0};
    std::size_t bits_used{0};
};

} // namespace adhoc

#endif // BRNG_TOOLS_SEED_SEQ_FILLER
