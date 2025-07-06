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

#ifndef BRNG_CHACHA_ENGINE
#define BRNG_CHACHA_ENGINE

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

template <class UIntType, std::size_t R, bool RFC7539 = false>
class chacha_engine final
    : public common_engine<UIntType, 32, chacha_engine<UIntType, R, RFC7539>> {

    static constexpr std::size_t w = 32;

    inline void increase_counter() {
        if constexpr (RFC7539) {
            ++this->state[12];
        } else {
            if (!++this->state[12]) {
                ++this->state[13];
            }
        }
    }

    inline void decrease_counter() {
        if constexpr (RFC7539) {
            --this->state[12];
        } else {
            if (!this->state[12]--) {
                --this->state[13];
            }
        }
    }

    inline void generate_block() {
        std::copy(this->state.begin(), this->state.end(),
                  this->output_cache.begin());

        auto quarterround =
            [&x = this->output_cache]<std::size_t a, std::size_t b,
                                      std::size_t c, std::size_t d>() {
                x[a] = x[a] + x[b];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[a] &= mask_result;
                }
                x[d] ^= x[a];
                x[d] = rotl<w>(x[d], 16);
                x[c] = x[c] + x[d];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[c] &= mask_result;
                }
                x[b] ^= x[c];
                x[b] = rotl<w>(x[b], 12);
                x[a] = x[a] + x[b];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[a] &= mask_result;
                }
                x[d] ^= x[a];
                x[d] = rotl<w>(x[d], 8);
                x[c] = x[c] + x[d];
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    x[c] &= mask_result;
                }
                x[b] ^= x[c];
                x[b] = rotl<w>(x[b], 7);
            };

        for (std::size_t round = 1; round < R; round += 2) {
            quarterround.template operator()<0, 4, 8, 12>();
            quarterround.template operator()<1, 5, 9, 13>();
            quarterround.template operator()<2, 6, 10, 14>();
            quarterround.template operator()<3, 7, 11, 15>();
            quarterround.template operator()<0, 5, 10, 15>();
            quarterround.template operator()<1, 6, 11, 12>();
            quarterround.template operator()<2, 7, 8, 13>();
            quarterround.template operator()<3, 4, 9, 14>();
        }
        if constexpr (static_cast<bool>(R & 1U)) {
            quarterround.template operator()<0, 5, 10, 15>();
            quarterround.template operator()<1, 6, 11, 12>();
            quarterround.template operator()<2, 7, 8, 13>();
            quarterround.template operator()<3, 4, 9, 14>();
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
    static constexpr UIntType mask_result = chacha_engine::max();

    chacha_engine() : chacha_engine(default_seed_1) {}

    chacha_engine(std::array<std::uint32_t, 8> const &key_input,
                  std::array<std::uint32_t, RFC7539 ? 3 : 2> const
                      &nonce_input = std::array<std::uint32_t, 2>{}) {
        constexpr std::string_view phrase{"expand 32-byte k"};
        constexpr auto constants = constants_from_phrase<phrase.size()>(phrase);

        std::copy(constants.begin(), constants.end(), state.begin());
        std::copy(key_input.begin(), key_input.end(),
                  state.begin() + constants.size());
        std::copy(nonce_input.begin(), nonce_input.end(),
                  state.begin() + constants.size() + key_input.size() +
                      (RFC7539 ? 1 : 2));

        this->generate_block();
    }

    chacha_engine(std::uint64_t seedval)
        : chacha_engine(std::array<std::uint32_t, 8>{seedval & mask_result,
                                                     (seedval >> 32U)}) {}

    void
    set_nonce(std::array<std::uint32_t, RFC7539 ? 3 : 2> const &nonce_input) {
        if constexpr (RFC7539) {
            this->state[13] = nonce_input[0];
            this->state[14] = nonce_input[1];
            this->state[15] = nonce_input[2];

        } else {
            this->state[14] = nonce_input[0];
            this->state[15] = nonce_input[1];
        }
    }

    inline auto operator*() const -> value_type {
        return this->output_cache[this->idx];
    };

    inline auto operator++() -> chacha_engine & {
        constexpr auto mask4 = mask<std::size_t>(4);

        ++this->idx;
        this->idx &= mask4;
        if (this->idx == 0) {
            this->increase_counter();
            this->generate_block();
        }

        return *this;
    }

    inline auto operator--() -> chacha_engine & {
        constexpr auto mask4 = mask<std::size_t>(4);

        if (this->idx == 0) {
            this->decrease_counter();
            this->generate_block();
        }

        --this->idx;
        this->idx &= mask4;

        return *this;
    }

    using common_engine<UIntType, w, chacha_engine>::operator++;
    using common_engine<UIntType, w, chacha_engine>::operator--;
    using common_engine<UIntType, w, chacha_engine>::operator==;

    auto operator==(const chacha_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->idx == rhs.idx) &&
               (this->output_cache == rhs.output_cache);
    }

  private:
    std::array<UIntType, 16> state{};
    std::size_t idx{0};
    std::array<UIntType, 16> output_cache{};
};

using chacha20 = chacha_engine<std::uint_fast32_t, 20>;

} // namespace adhoc

#endif // BRNG_CHACHA_ENGINE
