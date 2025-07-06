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

#ifndef BRNG_SALSA_ENGINE
#define BRNG_SALSA_ENGINE

#include "tools/bit.hpp"
#include "tools/common_engine.hpp"
#include "tools/mask.hpp"
#include "tools/string_utils.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>

namespace adhoc {

template <class UIntType, std::size_t R>
class salsa_engine final
    : public common_engine<UIntType, 32, salsa_engine<UIntType, R>> {

    static constexpr std::size_t w = 32;

    inline void increase_counter() {
        if (!++this->state[8]) {
            ++this->state[9];
        }
    }

    inline void decrease_counter() {
        if (!this->state[8]--) {
            --this->state[9];
        }
    }

    inline void generate_block() {
        std::copy(this->state.begin(), this->state.end(),
                  this->output_cache.begin());

        auto quarterround =
            [&x = this->output_cache]<std::size_t a, std::size_t b,
                                      std::size_t c, std::size_t d>() {
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[b] ^= rotl<w>((x[a] + x[d]) & mask_result, 7);
                    x[c] ^= rotl<w>((x[b] + x[a]) & mask_result, 9);
                    x[d] ^= rotl<w>((x[c] + x[b]) & mask_result, 13);
                    x[a] ^= rotl<w>((x[d] + x[c]) & mask_result, 18);
                } else {
                    x[b] ^= rotl<w>(x[a] + x[d], 7);
                    x[c] ^= rotl<w>(x[b] + x[a], 9);
                    x[d] ^= rotl<w>(x[c] + x[b], 13);
                    x[a] ^= rotl<w>(x[d] + x[c], 18);
                }
            };

        for (std::size_t round = 1; round < R; round += 2) {
            quarterround.template operator()<0, 4, 8, 12>();
            quarterround.template operator()<5, 9, 13, 1>();
            quarterround.template operator()<10, 14, 2, 6>();
            quarterround.template operator()<15, 3, 7, 11>();
            quarterround.template operator()<0, 1, 2, 3>();
            quarterround.template operator()<5, 6, 7, 4>();
            quarterround.template operator()<10, 11, 8, 9>();
            quarterround.template operator()<15, 12, 13, 14>();
        }
        if constexpr (static_cast<bool>(R & 1U)) {
            quarterround.template operator()<0, 4, 8, 12>();
            quarterround.template operator()<5, 9, 13, 1>();
            quarterround.template operator()<10, 14, 2, 6>();
            quarterround.template operator()<15, 3, 7, 11>();
        }
        for (std::size_t i = 0; i < 16; ++i) {
            this->output_cache[i] += this->state[i];
            if constexpr (w < std::numeric_limits<UIntType>::digits) {
                this->output_cache[i] &= mask_result;
            }
        }
    }

  public:
    using value_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::uint64_t default_seed_1 = 0xb504f333f9de6484UL;
    static constexpr UIntType mask_result = salsa_engine::max();

    salsa_engine() : salsa_engine(default_seed_1) {}

    salsa_engine(std::uint64_t seedval)
        : salsa_engine(std::array<std::uint32_t, 8>{seedval & mask_result,
                                                    (seedval >> 32U)}) {}

    salsa_engine(std::array<std::uint32_t, 8> const &key_input,
                 std::array<std::uint32_t, 2> const &nonce_input =
                     std::array<std::uint32_t, 2>{}) {

        constexpr std::string_view phrase{"expand 32-byte k"};
        constexpr auto constants = constants_from_phrase<phrase.size()>(phrase);

        this->state[0] = constants[0];
        this->state[5] = constants[1];
        this->state[10] = constants[2];
        this->state[15] = constants[3];
        for (int i = 0; i < 4; i++) {
            this->state[1 + i] = key_input[i];
        }
        for (int i = 0; i < 4; i++) {
            this->state[11 + i] = key_input[i + 4];
        }
        this->state[8] = 0;
        this->state[9] = 0;
        this->state[6] = nonce_input[0];
        this->state[7] = nonce_input[1];
        this->generate_block();
    }

    void set_nonce(std::array<std::uint32_t, 2> const &nonce_input) {
        this->state[6] = nonce_input[0];
        this->state[7] = nonce_input[1];
    }

    inline auto operator*() const -> value_type {
        return this->output_cache[this->idx];
    };

    inline auto operator++() -> salsa_engine & {
        constexpr auto mask4 = mask<std::size_t>(4);

        ++this->idx;
        this->idx &= mask4;
        if (this->idx == 0) {
            this->increase_counter();
            this->generate_block();
        }

        return *this;
    }

    inline auto operator--() -> salsa_engine & {
        constexpr auto mask4 = mask<std::size_t>(4);

        if (this->idx == 0) {
            this->decrease_counter();
            this->generate_block();
        }

        --this->idx;
        this->idx &= mask4;

        return *this;
    }

    using common_engine<UIntType, w, salsa_engine>::operator++;
    using common_engine<UIntType, w, salsa_engine>::operator--;
    using common_engine<UIntType, w, salsa_engine>::operator==;

    auto operator==(const salsa_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->idx == rhs.idx) &&
               (this->output_cache == rhs.output_cache);
    }

  private:
    std::array<UIntType, 16> state{};
    std::size_t idx{0};
    std::array<UIntType, 16> output_cache{};
};

using salsa20 = salsa_engine<std::uint_fast32_t, 20>;

} // namespace adhoc

#endif // BRNG_SALSA_ENGINE
