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

#ifndef BRNG_SFMT_ENGINE
#define BRNG_SFMT_ENGINE

#include "tools/circular_buffer.hpp"
#include "tools/mask.hpp"
#include "tools/uint128.hpp"
#include "tools/unshift.hpp"
// #include "tools/seed_seq_filler.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>

#include <cstring>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t SFMT_MEXP,
          std::size_t SFMT_POS1, std::size_t SFMT_SL1, std::size_t SFMT_SL2,
          std::size_t SFMT_SR1, std::size_t SFMT_SR2, std::uint32_t SFMT_MSK1,
          std::uint32_t SFMT_MSK2, std::uint32_t SFMT_MSK3,
          std::uint32_t SFMT_MSK4, std::uint32_t SFMT_PARITY1,
          std::uint32_t SFMT_PARITY2, std::uint32_t SFMT_PARITY3,
          std::uint32_t SFMT_PARITY4, bool use128_t = true>
class sfmt_engine final {
  private:
    static constexpr std::size_t SFMT_N = ((SFMT_MEXP + 127) / 128);
    static constexpr std::size_t SFMT_N32 = (SFMT_N * 4);
    static constexpr std::size_t SFMT_N64 = (SFMT_N * 2);

    using w128_t = std::array<std::uint32_t, 4>;

    static_assert(w <= std::numeric_limits<UIntType>::digits);

    void period_certification() {
        w128_t &first = this->state.at();
        constexpr std::array<std::uint32_t, 4> parity = {
            SFMT_PARITY1, SFMT_PARITY2, SFMT_PARITY3, SFMT_PARITY4};

        std::uint32_t inner = 0;
        for (std::size_t i = 0; i < 4; i++) {
            inner ^= first[i] & parity[i];
        }
        for (std::size_t i = 16; i > 0; i >>= 1) {
            inner ^= inner >> i;
        }
        inner &= 1U;
        // check OK
        if (inner == 1U) {
            return;
        }

        // check NG, and modification
        for (std::size_t i = 0; i < 4; i++) {
            std::uint32_t work = 1;
            for (std::size_t j = 0; j < 32; j++) {
                if ((work & parity[i]) != 0) {
                    first[i] ^= work;
                    return;
                }
                work <<= 1U;
            }
        }
    }

    void init_consistent() {
        this->period_certification();
        --this->state;
        this->operator++();
    }

    auto rshift128(w128_t const &in) {
        if constexpr (use128_t) {
            uint128 value = (static_cast<uint128>(in[3]) << 96U) |
                            (static_cast<uint128>(in[2]) << 64U) |
                            (static_cast<uint128>(in[1]) << 32U) |
                            static_cast<uint128>(in[0]);

            value >>= SFMT_SR2;

            return w128_t{static_cast<std::uint32_t>(value),
                          static_cast<std::uint32_t>(value >> 32U),
                          static_cast<std::uint32_t>(value >> 64U),
                          static_cast<std::uint32_t>(value >> 96U)};
        } else {
            std::uint64_t const in_high =
                (static_cast<std::uint64_t>(in[3]) << 32U) |
                static_cast<std::uint64_t>(in[2]);
            std::uint64_t const in_low =
                (static_cast<std::uint64_t>(in[1]) << 32U) |
                static_cast<std::uint64_t>(in[0]);

            std::uint64_t const out_high = in_high >> SFMT_SR2;
            std::uint64_t const out_low =
                (in_high << (64U - SFMT_SR2)) | (in_low >> SFMT_SR2);
            return w128_t{static_cast<std::uint32_t>(out_low),
                          static_cast<std::uint32_t>(out_low >> 32U),
                          static_cast<std::uint32_t>(out_high),
                          static_cast<std::uint32_t>(out_high >> 32U)};
        }
    }

    auto lshift128(w128_t const &in) {
        if constexpr (use128_t) {
            uint128 value = (static_cast<uint128>(in[3]) << 96U) |
                            (static_cast<uint128>(in[2]) << 64U) |
                            (static_cast<uint128>(in[1]) << 32U) |
                            static_cast<uint128>(in[0]);

            value <<= SFMT_SL2;

            return w128_t{static_cast<std::uint32_t>(value),
                          static_cast<std::uint32_t>(value >> 32U),
                          static_cast<std::uint32_t>(value >> 64U),
                          static_cast<std::uint32_t>(value >> 96U)};
        } else {
            std::uint64_t const in_high =
                (static_cast<std::uint64_t>(in[3]) << 32U) |
                static_cast<std::uint64_t>(in[2]);
            std::uint64_t const in_low =
                (static_cast<std::uint64_t>(in[1]) << 32U) |
                static_cast<std::uint64_t>(in[0]);

            std::uint64_t const out_high =
                (in_high << SFMT_SL2) | (in_low >> (64U - SFMT_SL2));
            std::uint64_t const out_low = in_low << SFMT_SL2;
            return w128_t{static_cast<std::uint32_t>(out_low),
                          static_cast<std::uint32_t>(out_low >> 32U),
                          static_cast<std::uint32_t>(out_high),
                          static_cast<std::uint32_t>(out_high >> 32U)};
        }
    }

  public:
    using difference_type = std::ptrdiff_t;
    using value_type = UIntType;

    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t state_size = SFMT_N32;
    static constexpr UIntType default_seed = 5489U;

    sfmt_engine() : sfmt_engine(default_seed) {}

    explicit sfmt_engine(std::uint32_t seed) {

        std::size_t i = 0;
        for (auto &arr : this->state.data()) {
            for (auto &val : arr) {
                val = seed;
                seed = 1812433253UL * (seed ^ (seed >> 30U)) + ++i;
            }
        }
        this->init_consistent();
    }

    template <std::size_t key_length>
    explicit sfmt_engine(std::array<std::uint32_t, key_length> init_key) {

        std::uint32_t rr;
        std::uint32_t *psfmt32 = &this->state.data()[0][0];

        constexpr std::size_t size = SFMT_N * 4;
        constexpr std::size_t lag = size >= 623  ? 11
                                    : size >= 68 ? 7
                                    : size >= 39 ? 5
                                                 : 3;
        constexpr std::size_t mid = (size - lag) / 2;

        for (auto &val : this->state.data()) {
            for (auto &val2 : val) {
                val2 = 2341178251;
            }
        }

        constexpr std::size_t count = std::max(SFMT_N32, key_length + 1U);

        for (std::size_t j = 0; j < count; ++j) {
            std::size_t const i = j % SFMT_N32;
            std::uint32_t const x = psfmt32[i] ^ psfmt32[(i + mid) % SFMT_N32] ^
                                    psfmt32[(i + SFMT_N32 - 1) % SFMT_N32];

            rr = (x ^ (x >> 27U)) * static_cast<std::uint32_t>(1664525UL);
            psfmt32[(i + mid) % SFMT_N32] += rr;

            if (j == 0) {
                rr += key_length;
            } else {
                rr += i;
                if ((j - 1) < key_length) {
                    rr += init_key[j - 1];
                }
            }

            psfmt32[(i + mid + lag) % SFMT_N32] += rr;
            psfmt32[i] = rr;
        }

        for (std::size_t i = 0; i < SFMT_N32; i++) {
            std::uint32_t const x = psfmt32[i] + psfmt32[(i + mid) % SFMT_N32] +
                                    psfmt32[(i + SFMT_N32 - 1) % SFMT_N32];
            rr = (x ^ (x >> 27U)) * static_cast<std::uint32_t>(1566083941UL);

            psfmt32[(i + mid) % SFMT_N32] ^= rr;
            rr -= i;
            psfmt32[(i + mid + lag) % SFMT_N32] ^= rr;
            psfmt32[i] = rr;
        }
        this->init_consistent();
    }

    template <class SeedSeq> explicit sfmt_engine(SeedSeq &seq) {

        std::array<std::uint32_t, SFMT_N32> in{};
        seq.generate(in.begin(), in.end());
        std::size_t in_idx = 0;
        for (auto &arr : this->state.data()) {
            for (auto &val : arr) {
                val = in[in_idx++];
            }
        }

        this->init_consistent();
    }

    inline auto operator*() const -> value_type {
        return this->state.at()[this->idx];
    };

    inline auto operator++() -> sfmt_engine & {
        ++this->idx;
        if (this->idx == 4) {
            ++this->state;
            this->idx = 0;

            w128_t &current = this->state.at();
            w128_t const &right_shifted = this->state.template at<SFMT_POS1>();
            w128_t const &left_shifted = this->state.template at<SFMT_N - 1>();
            w128_t const &right_shifted_128 =
                this->state.template at<SFMT_N - 2>();
            auto const x = lshift128(current);
            auto const y = rshift128(right_shifted_128);
            current[0] ^= x[0] ^ ((right_shifted[0] >> SFMT_SR1) & SFMT_MSK1) ^
                          y[0] ^ (left_shifted[0] << SFMT_SL1);
            current[1] ^= x[1] ^ ((right_shifted[1] >> SFMT_SR1) & SFMT_MSK2) ^
                          y[1] ^ (left_shifted[1] << SFMT_SL1);
            current[2] ^= x[2] ^ ((right_shifted[2] >> SFMT_SR1) & SFMT_MSK3) ^
                          y[2] ^ (left_shifted[2] << SFMT_SL1);
            current[3] ^= x[3] ^ ((right_shifted[3] >> SFMT_SR1) & SFMT_MSK4) ^
                          y[3] ^ (left_shifted[3] << SFMT_SL1);
        }
        return *this;
    }

    inline auto operator--() -> sfmt_engine & {
        if (this->idx == 0) {
            w128_t &current = this->state.at();
            w128_t const &right_shifted = this->state.template at<SFMT_POS1>();
            w128_t const &left_shifted = this->state.template at<SFMT_N - 1>();
            w128_t const &right_shifted_128 =
                this->state.template at<SFMT_N - 2>();
            auto const y = rshift128(right_shifted_128);

            current[0] ^= ((right_shifted[0] >> SFMT_SR1) & SFMT_MSK1) ^ y[0] ^
                          (left_shifted[0] << SFMT_SL1);
            current[1] ^= ((right_shifted[1] >> SFMT_SR1) & SFMT_MSK2) ^ y[1] ^
                          (left_shifted[1] << SFMT_SL1);
            current[2] ^= ((right_shifted[2] >> SFMT_SR1) & SFMT_MSK3) ^ y[2] ^
                          (left_shifted[2] << SFMT_SL1);
            current[3] ^= ((right_shifted[3] >> SFMT_SR1) & SFMT_MSK4) ^ y[3] ^
                          (left_shifted[3] << SFMT_SL1);

            uint128 value = (static_cast<uint128>(current[3]) << 96U) |
                            (static_cast<uint128>(current[2]) << 64U) |
                            (static_cast<uint128>(current[1]) << 32U) |
                            static_cast<uint128>(current[0]);

            value = unshift_left_xor<128, SFMT_SL2>(value);

            current[0] = static_cast<std::uint32_t>(value);
            current[1] = static_cast<std::uint32_t>(value >> 32U);
            current[2] = static_cast<std::uint32_t>(value >> 64U);
            current[3] = static_cast<std::uint32_t>(value >> 96U);

            this->idx = 4;
            --this->state;
        }
        --this->idx;
        return *this;
    }

    // creates a copy. should be avoided!
    inline auto operator++(int) -> sfmt_engine {
        auto const tmp = *this;
        ++*this;
        return tmp;
    }

    // creates a copy. should be avoided!
    inline auto operator--(int) -> sfmt_engine {
        auto const tmp = *this;
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
    static constexpr auto max() -> UIntType { return mask<UIntType>(w); }

    auto operator==(const sfmt_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->idx == rhs.idx);
    }

    auto operator!=(const sfmt_engine &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

    constexpr auto operator==(std::default_sentinel_t /*unused*/) const
        -> bool {
        return false;
    }

  private:
    circular_buffer<w128_t, SFMT_N> state{};
    std::size_t idx{4 - 1};
};

using sfmt607 = sfmt_engine<std::uint_fast32_t, 32, 607, 2, 15, 24U, 13, 24U,
                            0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                            0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>;

using sfmt1279 =
    sfmt_engine<std::uint_fast32_t, 32, 1279, 7, 14, 24U, 5, 8U, 0xf7fefffdU,
                0x7fefcfffU, 0xaff3ef3fU, 0xb5ffff7fU, 0x00000001U, 0x00000000U,
                0x00000000U, 0x20000000U>;

using sfmt2281 =
    sfmt_engine<std::uint_fast32_t, 32, 2281, 12, 19, 8U, 5, 8U, 0xbff7ffbfU,
                0xfdfffffeU, 0xf7ffef7fU, 0xf2f7cbbfU, 0x00000001U, 0x00000000U,
                0x00000000U, 0x41dfa600U>;

using sfmt4253 =
    sfmt_engine<std::uint_fast32_t, 32, 4253, 17, 20, 8U, 7, 8U, 0x9f7bffffU,
                0x9fffff5fU, 0x3efffffbU, 0xfffff7bbU, 0xa8000001U, 0xaf5390a3U,
                0xb740b3f8U, 0x6c11486dU>;

using sfmt11213 =
    sfmt_engine<std::uint_fast32_t, 32, 11213, 68, 14, 24U, 7, 24U, 0xeffff7fbU,
                0xffffffefU, 0xdfdfbfffU, 0x7fffdbfdU, 0x00000001U, 0x00000000U,
                0xe8148000U, 0xd0c7afa3U>;

using sfmt19937 =
    sfmt_engine<std::uint_fast32_t, 32, 19937, 122, 18, 8U, 11, 8U, 0xdfffffefU,
                0xddfecb7fU, 0xbffaffffU, 0xbffffff6U, 0x00000001U, 0x00000000U,
                0x00000000U, 0x13c9e684U>;

using sfmt44497 =
    sfmt_engine<std::uint_fast32_t, 32, 44497, 330, 5, 24U, 9, 24U, 0xeffffffbU,
                0xdfbebfffU, 0xbfbf7befU, 0x9ffd7bffU, 0x00000001U, 0x00000000U,
                0xa3ac4000U, 0xecc1327aU>;

using sfmt86243 =
    sfmt_engine<std::uint_fast32_t, 32, 86243, 366, 6, 56U, 19, 8U, 0xfdbffbffU,
                0xbff7ff3fU, 0xfd77efffU, 0xbf9ff3ffU, 0x00000001U, 0x00000000U,
                0x00000000U, 0xe9528d85U>;

using sfmt132049 =
    sfmt_engine<std::uint_fast32_t, 32, 132049, 110, 19, 8U, 21, 8U,
                0xffffbb5fU, 0xfb6ebf95U, 0xfffefffaU, 0xcff77fffU, 0x00000001U,
                0x00000000U, 0xcb520000U, 0xc7e91c7dU>;

using sfmt216091 =
    sfmt_engine<std::uint_fast32_t, 32, 216091, 627, 11, 24U, 10, 8U,
                0xbff7bff7U, 0xbfffffffU, 0xbffffa7fU, 0xffddfbfbU, 0xf8000001U,
                0x89e80709U, 0x3bd2b64bU, 0x0c64b1e4U>;

} // namespace adhoc

#endif // BRNG_SFMT_ENGINE
