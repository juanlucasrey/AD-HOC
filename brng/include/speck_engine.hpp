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

#ifndef BRNG_SPECK_ENGINE
#define BRNG_SPECK_ENGINE

#include "tools/common_engine.hpp"
#include "tools/mask.hpp"
#include "tools/seed_seq_filler.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace adhoc {

template <class UIntType, std::size_t R>
class speck_engine final
    : public common_engine<UIntType, 64, speck_engine<UIntType, R>> {

    static constexpr std::size_t w = 64;

    inline void increase_counter() {
        if (!++this->ctr[0]) {
            ++this->ctr[1];
        }
    }

    inline void decrease_counter() {
        if (!this->ctr[0]--) {
            --this->ctr[1];
        }
    }

    inline void generate_block() {
        this->output_cache[0] = this->ctr[0];
        this->output_cache[1] = this->ctr[1];
        for (std::size_t i = 0; i < R; i++) {
            this->output_cache[1] =
                (std::rotr(this->output_cache[1], 8) + this->output_cache[0]) ^
                (this->state[i]);
            this->output_cache[0] =
                std::rotl(this->output_cache[0], 3) ^ this->output_cache[1];
        }
    }

  public:
    using value_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::uint64_t default_seed_1 = 0xb504f333f9de6484UL;
    static constexpr UIntType mask_result = speck_engine::max();

    speck_engine() : speck_engine(default_seed_1) {}

    speck_engine(std::uint64_t input1, std::uint64_t input2 = 0,
                 std::uint64_t input3 = 0, std::uint64_t input4 = 0) {

        std::uint64_t i;
        std::uint64_t D = input4;
        std::uint64_t C = input3;
        std::uint64_t B = input2;
        std::uint64_t A = input1;
        for (std::size_t i = 0; i < 33; i += 3) {
            this->state[i] = A;
            B = (std::rotr(B, 8) + A) ^ (i);
            A = std::rotl(A, 3) ^ B;
            this->state[i + 1] = A;
            C = (std::rotr(C, 8) + A) ^ (i + 1);
            A = std::rotl(A, 3) ^ C;
            this->state[i + 2] = A;
            D = (std::rotr(D, 8) + A) ^ (i + 2);
            A = std::rotl(A, 3) ^ D;
        }
        this->state[33] = A;
        this->generate_block();
    }

    template <class SeedSeq> explicit speck_engine(SeedSeq &seq) {
        seed_seq_filler<w, 34> seq_filler(seq);
        seq_filler.template generate<UIntType>(this->state.begin(),
                                               this->state.end());
        this->generate_block();
    }

    inline auto operator*() const -> value_type {
        return this->output_cache[this->idx];
    };

    inline auto operator++() -> speck_engine & {
        constexpr auto mask1 = mask<std::size_t>(1);

        ++this->idx;
        this->idx &= mask1;
        if (this->idx == 0) {
            this->increase_counter();
            this->generate_block();
        }

        return *this;
    }

    inline auto operator--() -> speck_engine & {
        constexpr auto mask1 = mask<std::size_t>(1);

        if (this->idx == 0) {
            this->decrease_counter();
            this->generate_block();
        }

        --this->idx;
        this->idx &= mask1;

        return *this;
    }

    using common_engine<UIntType, w, speck_engine>::operator++;
    using common_engine<UIntType, w, speck_engine>::operator--;
    using common_engine<UIntType, w, speck_engine>::operator==;

    auto operator==(const speck_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->ctr == rhs.ctr) &&
               (this->output_cache == rhs.output_cache) &&
               (this->idx == rhs.idx);
    }

  private:
    std::array<std::uint64_t, 34> state{};
    std::array<std::uint64_t, 2> ctr{};
    std::array<std::uint64_t, 2> output_cache{};
    std::size_t idx{0};
};

using speck34 = speck_engine<std::uint64_t, 34>;

} // namespace adhoc

#endif // BRNG_SPECK_ENGINE
