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
#include "tools/common_engine.hpp"
#include "tools/mask.hpp"
#include "tools/seed_seq_filler.hpp"
#include "tools/uint128.hpp"
#include "tools/unshift.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t mexp, std::size_t pos1,
          std::size_t sl1, std::size_t sl2, std::size_t sr1, std::size_t sr2,
          std::uint32_t msk1, std::uint32_t msk2, std::uint32_t msk3,
          std::uint32_t msk4, std::uint32_t parity1, std::uint32_t parity2,
          std::uint32_t parity3, std::uint32_t parity4, bool use128_t = true>
class sfmt_engine final
    : public common_engine<UIntType, w,
                           sfmt_engine<UIntType, w, mexp, pos1, sl1, sl2, sr1,
                                       sr2, msk1, msk2, msk3, msk4, parity1,
                                       parity2, parity3, parity4, use128_t>> {
  private:
    static_assert(w == 32 || w == 64);
    static constexpr std::size_t n = ((mexp + 127) / 128);
    static constexpr std::size_t sub_size = 128 / w;
    static constexpr UIntType mask_result = sfmt_engine::max();

    using w128_t = std::array<UIntType, sub_size>;

    auto get_val32(std::size_t i) -> std::uint32_t {
        auto &state_data = this->state.data();
        std::size_t const state_idx = i / 4;
        std::size_t const sub_idx = i % 4;
        if constexpr (w == 32) {
            return state_data[state_idx][sub_idx];
        } else {
            std::size_t const sub_idx_64 = sub_idx / 2;
            std::size_t const sub_idx_64_half = sub_idx % 2;
            UIntType val = state_data[state_idx][sub_idx_64];
            if (sub_idx_64_half == 1) {
                val >>= 32U;
            }
            return static_cast<std::uint32_t>(val);
        }
    };

    void set_val32(std::size_t i, std::uint32_t r) {
        auto &state_data = this->state.data();
        std::size_t const state_idx = i / 4;
        std::size_t const sub_idx = i % 4;
        if constexpr (w == 32) {
            state_data[state_idx][sub_idx] = r;
        } else {
            std::size_t const sub_idx_64 = sub_idx / 2;
            std::size_t const sub_idx_64_half = sub_idx % 2;
            UIntType &val = state_data[state_idx][sub_idx_64];
            auto r64 = static_cast<UIntType>(r);
            if (sub_idx_64_half == 1) {
                constexpr auto mask_lower_half = mask<UIntType>(32);
                r64 <<= 32U;
                val = r64 | (val & mask_lower_half);
            } else {
                constexpr auto mask_upper_half = mask<UIntType>(32, 32);
                val = r64 | (val & mask_upper_half);
            }
        }
    };

    void period_certification() {
        w128_t &first = this->state.at();
        constexpr std::array<std::uint32_t, 4> parity = {parity1, parity2,
                                                         parity3, parity4};

        std::uint32_t inner = 0;
        for (std::size_t i = 0; i < 4; i++) {
            inner ^= get_val32(i) & parity[i];
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
                    set_val32(i, get_val32(i) ^ work);
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

    inline auto rshift128(w128_t const &in) {
        if constexpr (use128_t) {
            if constexpr (w == 32) {
                uint128 value = (static_cast<uint128>(in[3]) << 96U) |
                                (static_cast<uint128>(in[2]) << 64U) |
                                (static_cast<uint128>(in[1]) << 32U) |
                                static_cast<uint128>(in[0]);

                value >>= sr2;

                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{
                        static_cast<UIntType>(value) & mask_result,
                        static_cast<UIntType>(value >> 32U) & mask_result,
                        static_cast<UIntType>(value >> 64U) & mask_result,
                        static_cast<UIntType>(value >> 96U) & mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(value),
                                  static_cast<UIntType>(value >> 32U),
                                  static_cast<UIntType>(value >> 64U),
                                  static_cast<UIntType>(value >> 96U)};
                }
            } else {
                uint128 value = (static_cast<uint128>(in[1]) << 64U) |
                                static_cast<uint128>(in[0]);

                value >>= sr2;

                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{static_cast<UIntType>(value) & mask_result,
                                  static_cast<UIntType>(value >> 64U) &
                                      mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(value),
                                  static_cast<UIntType>(value >> 64U)};
                }
            }
        } else {
            if constexpr (w == 32) {
                std::uint64_t const in_high =
                    (static_cast<std::uint64_t>(in[3]) << 32U) |
                    static_cast<std::uint64_t>(in[2]);
                std::uint64_t const in_low =
                    (static_cast<std::uint64_t>(in[1]) << 32U) |
                    static_cast<std::uint64_t>(in[0]);

                std::uint64_t const out_high = in_high >> sr2;
                std::uint64_t const out_low =
                    (in_high << (64U - sr2)) | (in_low >> sr2);
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{
                        static_cast<UIntType>(out_low) & mask_result,
                        static_cast<UIntType>(out_low >> 32U) & mask_result,
                        static_cast<UIntType>(out_high) & mask_result,
                        static_cast<UIntType>(out_high >> 32U) & mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(out_low),
                                  static_cast<UIntType>(out_low >> 32U),
                                  static_cast<UIntType>(out_high),
                                  static_cast<UIntType>(out_high >> 32U)};
                }
            } else {
                UIntType const out_high = in[1] >> sr2;
                UIntType const out_low =
                    (in[1] << (64U - sr2)) | (in[0] >> sr2);
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{static_cast<UIntType>(out_low) & mask_result,
                                  static_cast<UIntType>(out_high) &
                                      mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(out_low),
                                  static_cast<UIntType>(out_high)};
                }
            }
        }
    }

    inline auto lshift128(w128_t const &in) {
        if constexpr (use128_t) {
            if constexpr (w == 32) {
                uint128 value = (static_cast<uint128>(in[3]) << 96U) |
                                (static_cast<uint128>(in[2]) << 64U) |
                                (static_cast<uint128>(in[1]) << 32U) |
                                static_cast<uint128>(in[0]);

                value <<= sl2;

                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{
                        static_cast<UIntType>(value) & mask_result,
                        static_cast<UIntType>(value >> 32U) & mask_result,
                        static_cast<UIntType>(value >> 64U) & mask_result,
                        static_cast<UIntType>(value >> 96U) & mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(value),
                                  static_cast<UIntType>(value >> 32U),
                                  static_cast<UIntType>(value >> 64U),
                                  static_cast<UIntType>(value >> 96U)};
                }
            } else {
                uint128 value = (static_cast<uint128>(in[1]) << 64U) |
                                static_cast<uint128>(in[0]);

                value <<= sl2;

                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{static_cast<UIntType>(value) & mask_result,
                                  static_cast<UIntType>(value >> 64U) &
                                      mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(value),
                                  static_cast<UIntType>(value >> 64U)};
                }
            }
        } else {
            if constexpr (w == 32) {
                std::uint64_t const in_high =
                    (static_cast<std::uint64_t>(in[3]) << 32U) |
                    static_cast<std::uint64_t>(in[2]);
                std::uint64_t const in_low =
                    (static_cast<std::uint64_t>(in[1]) << 32U) |
                    static_cast<std::uint64_t>(in[0]);

                std::uint64_t const out_high =
                    (in_high << sl2) | (in_low >> (64U - sl2));
                std::uint64_t const out_low = in_low << sl2;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{
                        static_cast<UIntType>(out_low) & mask_result,
                        static_cast<UIntType>(out_low >> 32U) & mask_result,
                        static_cast<UIntType>(out_high) & mask_result,
                        static_cast<UIntType>(out_high >> 32U) & mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(out_low),
                                  static_cast<UIntType>(out_low >> 32U),
                                  static_cast<UIntType>(out_high),
                                  static_cast<UIntType>(out_high >> 32U)};
                }
            } else {
                UIntType const out_high =
                    (in[1] << sl2) | (in[0] >> (64U - sl2));
                UIntType const out_low = in[0] << sl2;
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    return w128_t{static_cast<UIntType>(out_low) & mask_result,
                                  static_cast<UIntType>(out_high) &
                                      mask_result};
                } else {
                    return w128_t{static_cast<UIntType>(out_low),
                                  static_cast<UIntType>(out_high)};
                }
            }
        }
    }

  public:
    using value_type = UIntType;
    static constexpr std::size_t word_size = w;
    static constexpr UIntType default_seed = 5489U;

    sfmt_engine() : sfmt_engine(default_seed) {}

    explicit sfmt_engine(std::uint32_t seed) {
        constexpr std::size_t size = n * 4;
        for (std::size_t i = 0; i < size; ++i) {
            set_val32(i, seed);
            seed = 1812433253UL * (seed ^ (seed >> 30U)) + i + 1;
        }

        this->init_consistent();
    }

    template <std::size_t key_length>
    explicit sfmt_engine(std::array<std::uint32_t, key_length> init_key) {
        constexpr std::size_t size_32 = n * 4;
        constexpr std::size_t lag = size_32 >= 623  ? 11
                                    : size_32 >= 68 ? 7
                                    : size_32 >= 39 ? 5
                                                    : 3;
        constexpr std::size_t mid = (size_32 - lag) / 2;

        for (std::size_t i = 0; i < size_32; i++) {
            set_val32(i, 2341178251U);
        }

        constexpr std::size_t count = std::max(size_32, key_length + 1U);

        for (std::size_t j = 0; j < count; ++j) {
            std::size_t const i = j % size_32;
            std::uint32_t r = get_val32(i) ^ get_val32((i + mid) % size_32) ^
                              get_val32((i + size_32 - 1) % size_32);

            r = (r ^ (r >> 27U)) * 1664525UL;
            set_val32((i + mid) % size_32, get_val32((i + mid) % size_32) + r);

            if (j == 0) {
                r += key_length;
            } else {
                r += i;
                if ((j - 1) < key_length) {
                    r += init_key[j - 1];
                }
            }

            set_val32((i + mid + lag) % size_32,
                      get_val32((i + mid + lag) % size_32) + r);
            set_val32(i, r);
        }

        for (std::size_t i = 0; i < size_32; i++) {
            std::uint32_t r = get_val32(i) + get_val32((i + mid) % size_32) +
                              get_val32((i + size_32 - 1) % size_32);
            r = (r ^ (r >> 27U)) * 1566083941UL;

            set_val32((i + mid) % size_32, get_val32((i + mid) % size_32) ^ r);
            r -= i;
            set_val32((i + mid + lag) % size_32,
                      get_val32((i + mid + lag) % size_32) ^ r);
            set_val32(i, r);
        }
        this->init_consistent();
    }

    template <class SeedSeq> explicit sfmt_engine(SeedSeq &seq) {

        constexpr std::size_t size = n * sub_size;
        seed_seq_filler<w, size> seq_filler(seq);

        for (auto &arr : this->state.data()) {
            for (auto &val : arr) {
                val = seq_filler.template get<UIntType>();
            }
        }

        this->init_consistent();
    }

    inline auto operator*() const -> value_type {
        return this->state.at()[this->idx];
    };

    inline auto operator++() -> sfmt_engine & {
        ++this->idx;
        if (this->idx == sub_size) {
            ++this->state;
            this->idx = 0;

            auto &current = this->state.at();
            auto const &right_shifted = this->state.template at<pos1>();
            auto const &left_shifted = this->state.template at<n - 1>();
            auto const &right_shifted_128 = this->state.template at<n - 2>();
            auto const x = lshift128(current);
            auto const y = rshift128(right_shifted_128);
            if constexpr (w == 32) {
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    current[0] ^= x[0] ^ ((right_shifted[0] >> sr1) & msk1) ^
                                  y[0] ^
                                  ((left_shifted[0] << sl1) & mask_result);
                    current[1] ^= x[1] ^ ((right_shifted[1] >> sr1) & msk2) ^
                                  y[1] ^
                                  ((left_shifted[1] << sl1) & mask_result);
                    current[2] ^= x[2] ^ ((right_shifted[2] >> sr1) & msk3) ^
                                  y[2] ^
                                  ((left_shifted[2] << sl1) & mask_result);
                    current[3] ^= x[3] ^ ((right_shifted[3] >> sr1) & msk4) ^
                                  y[3] ^
                                  ((left_shifted[3] << sl1) & mask_result);
                } else {
                    current[0] ^= x[0] ^ ((right_shifted[0] >> sr1) & msk1) ^
                                  y[0] ^ (left_shifted[0] << sl1);
                    current[1] ^= x[1] ^ ((right_shifted[1] >> sr1) & msk2) ^
                                  y[1] ^ (left_shifted[1] << sl1);
                    current[2] ^= x[2] ^ ((right_shifted[2] >> sr1) & msk3) ^
                                  y[2] ^ (left_shifted[2] << sl1);
                    current[3] ^= x[3] ^ ((right_shifted[3] >> sr1) & msk4) ^
                                  y[3] ^ (left_shifted[3] << sl1);
                }
            } else {
                constexpr auto msk12 = (static_cast<UIntType>(msk2) << 32U) |
                                       (static_cast<UIntType>(msk1) &
                                        ~mask<UIntType>(sr1, 32 - sr1));
                constexpr auto msk34 = (static_cast<UIntType>(msk4) << 32U) |
                                       (static_cast<UIntType>(msk3) &
                                        ~mask<UIntType>(sr1, 32 - sr1));
                constexpr auto msksl = ~mask<UIntType>(sl1, 32) & mask_result;

                current[0] ^= x[0] ^ ((right_shifted[0] >> sr1) & msk12) ^
                              y[0] ^ ((left_shifted[0] << sl1) & msksl);
                current[1] ^= x[1] ^ ((right_shifted[1] >> sr1) & msk34) ^
                              y[1] ^ ((left_shifted[1] << sl1) & msksl);
            }
        }
        return *this;
    }

    inline auto operator--() -> sfmt_engine & {
        if (this->idx == 0) {
            auto &current = this->state.at();
            auto const &right_shifted = this->state.template at<pos1>();
            auto const &left_shifted = this->state.template at<n - 1>();
            auto const &right_shifted_128 = this->state.template at<n - 2>();
            auto const y = rshift128(right_shifted_128);

            if constexpr (w == 32) {
                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    current[0] ^= ((right_shifted[0] >> sr1) & msk1) ^ y[0] ^
                                  ((left_shifted[0] << sl1) & mask_result);
                    current[1] ^= ((right_shifted[1] >> sr1) & msk2) ^ y[1] ^
                                  ((left_shifted[1] << sl1) & mask_result);
                    current[2] ^= ((right_shifted[2] >> sr1) & msk3) ^ y[2] ^
                                  ((left_shifted[2] << sl1) & mask_result);
                    current[3] ^= ((right_shifted[3] >> sr1) & msk4) ^ y[3] ^
                                  ((left_shifted[3] << sl1) & mask_result);
                } else {
                    current[0] ^= ((right_shifted[0] >> sr1) & msk1) ^ y[0] ^
                                  (left_shifted[0] << sl1);
                    current[1] ^= ((right_shifted[1] >> sr1) & msk2) ^ y[1] ^
                                  (left_shifted[1] << sl1);
                    current[2] ^= ((right_shifted[2] >> sr1) & msk3) ^ y[2] ^
                                  (left_shifted[2] << sl1);
                    current[3] ^= ((right_shifted[3] >> sr1) & msk4) ^ y[3] ^
                                  (left_shifted[3] << sl1);
                }

                uint128 value = (static_cast<uint128>(current[3]) << 96U) |
                                (static_cast<uint128>(current[2]) << 64U) |
                                (static_cast<uint128>(current[1]) << 32U) |
                                static_cast<uint128>(current[0]);

                value = unshift_left_xor<128, sl2>(value);

                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    current[0] = static_cast<UIntType>(value) & mask_result;
                    current[1] =
                        static_cast<UIntType>(value >> 32U) & mask_result;
                    current[2] =
                        static_cast<UIntType>(value >> 64U) & mask_result;
                    current[3] =
                        static_cast<UIntType>(value >> 96U) & mask_result;
                } else {
                    current[0] = static_cast<UIntType>(value);
                    current[1] = static_cast<UIntType>(value >> 32U);
                    current[2] = static_cast<UIntType>(value >> 64U);
                    current[3] = static_cast<UIntType>(value >> 96U);
                }
            } else {
                constexpr auto msk12 = (static_cast<UIntType>(msk2) << 32U) |
                                       (static_cast<UIntType>(msk1) &
                                        ~mask<UIntType>(sr1, 32 - sr1));
                constexpr auto msk34 = (static_cast<UIntType>(msk4) << 32U) |
                                       (static_cast<UIntType>(msk3) &
                                        ~mask<UIntType>(sr1, 32 - sr1));
                constexpr auto msksl = ~mask<UIntType>(sl1, 32) & mask_result;
                current[0] ^= ((right_shifted[0] >> sr1) & msk12) ^ y[0] ^
                              ((left_shifted[0] << sl1) & msksl);
                current[1] ^= ((right_shifted[1] >> sr1) & msk34) ^ y[1] ^
                              ((left_shifted[1] << sl1) & msksl);

                uint128 value = (static_cast<uint128>(current[1]) << 64U) |
                                static_cast<uint128>(current[0]);

                value = unshift_left_xor<128, sl2>(value);

                if constexpr (w < std::numeric_limits<UIntType>::digits) {
                    current[0] = static_cast<UIntType>(value) & mask_result;
                    current[1] =
                        static_cast<UIntType>(value >> 64U) & mask_result;
                } else {
                    current[0] = static_cast<UIntType>(value);
                    current[1] = static_cast<UIntType>(value >> 64U);
                }
            }

            this->idx = sub_size;
            --this->state;
        }
        --this->idx;
        return *this;
    }

    using common_engine<UIntType, w, sfmt_engine>::operator++;
    using common_engine<UIntType, w, sfmt_engine>::operator--;

    auto operator==(const sfmt_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->idx == rhs.idx);
    }

  private:
    circular_buffer<w128_t, n> state{};
    std::size_t idx{sub_size - 1};
};

using sfmt607 = sfmt_engine<std::uint_fast32_t, 32, 607, 2, 15, 24U, 13, 24U,
                            0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                            0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>;

using sfmt607_64 =
    sfmt_engine<std::uint_fast64_t, 64, 607, 2, 15, 24U, 13, 24U, 0xfdff37ffU,
                0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU, 0x00000001U, 0x00000000U,
                0x00000000U, 0x5986f054U>;

using sfmt1279 =
    sfmt_engine<std::uint_fast32_t, 32, 1279, 7, 14, 24U, 5, 8U, 0xf7fefffdU,
                0x7fefcfffU, 0xaff3ef3fU, 0xb5ffff7fU, 0x00000001U, 0x00000000U,
                0x00000000U, 0x20000000U>;

using sfmt1279_64 =
    sfmt_engine<std::uint_fast64_t, 64, 1279, 7, 14, 24U, 5, 8U, 0xf7fefffdU,
                0x7fefcfffU, 0xaff3ef3fU, 0xb5ffff7fU, 0x00000001U, 0x00000000U,
                0x00000000U, 0x20000000U>;

using sfmt2281 =
    sfmt_engine<std::uint_fast32_t, 32, 2281, 12, 19, 8U, 5, 8U, 0xbff7ffbfU,
                0xfdfffffeU, 0xf7ffef7fU, 0xf2f7cbbfU, 0x00000001U, 0x00000000U,
                0x00000000U, 0x41dfa600U>;

using sfmt2281_64 =
    sfmt_engine<std::uint_fast64_t, 64, 2281, 12, 19, 8U, 5, 8U, 0xbff7ffbfU,
                0xfdfffffeU, 0xf7ffef7fU, 0xf2f7cbbfU, 0x00000001U, 0x00000000U,
                0x00000000U, 0x41dfa600U>;

using sfmt4253 =
    sfmt_engine<std::uint_fast32_t, 32, 4253, 17, 20, 8U, 7, 8U, 0x9f7bffffU,
                0x9fffff5fU, 0x3efffffbU, 0xfffff7bbU, 0xa8000001U, 0xaf5390a3U,
                0xb740b3f8U, 0x6c11486dU>;

using sfmt4253_64 =
    sfmt_engine<std::uint_fast64_t, 64, 4253, 17, 20, 8U, 7, 8U, 0x9f7bffffU,
                0x9fffff5fU, 0x3efffffbU, 0xfffff7bbU, 0xa8000001U, 0xaf5390a3U,
                0xb740b3f8U, 0x6c11486dU>;

using sfmt11213 =
    sfmt_engine<std::uint_fast32_t, 32, 11213, 68, 14, 24U, 7, 24U, 0xeffff7fbU,
                0xffffffefU, 0xdfdfbfffU, 0x7fffdbfdU, 0x00000001U, 0x00000000U,
                0xe8148000U, 0xd0c7afa3U>;

using sfmt11213_64 =
    sfmt_engine<std::uint_fast64_t, 64, 11213, 68, 14, 24U, 7, 24U, 0xeffff7fbU,
                0xffffffefU, 0xdfdfbfffU, 0x7fffdbfdU, 0x00000001U, 0x00000000U,
                0xe8148000U, 0xd0c7afa3U>;

using sfmt19937 =
    sfmt_engine<std::uint_fast32_t, 32, 19937, 122, 18, 8U, 11, 8U, 0xdfffffefU,
                0xddfecb7fU, 0xbffaffffU, 0xbffffff6U, 0x00000001U, 0x00000000U,
                0x00000000U, 0x13c9e684U>;

using sfmt19937_64 =
    sfmt_engine<std::uint_fast64_t, 64, 19937, 122, 18, 8U, 11, 8U, 0xdfffffefU,
                0xddfecb7fU, 0xbffaffffU, 0xbffffff6U, 0x00000001U, 0x00000000U,
                0x00000000U, 0x13c9e684U>;

using sfmt44497 =
    sfmt_engine<std::uint_fast32_t, 32, 44497, 330, 5, 24U, 9, 24U, 0xeffffffbU,
                0xdfbebfffU, 0xbfbf7befU, 0x9ffd7bffU, 0x00000001U, 0x00000000U,
                0xa3ac4000U, 0xecc1327aU>;

using sfmt44497_64 =
    sfmt_engine<std::uint_fast64_t, 64, 44497, 330, 5, 24U, 9, 24U, 0xeffffffbU,
                0xdfbebfffU, 0xbfbf7befU, 0x9ffd7bffU, 0x00000001U, 0x00000000U,
                0xa3ac4000U, 0xecc1327aU>;

using sfmt86243 =
    sfmt_engine<std::uint_fast32_t, 32, 86243, 366, 6, 56U, 19, 8U, 0xfdbffbffU,
                0xbff7ff3fU, 0xfd77efffU, 0xbf9ff3ffU, 0x00000001U, 0x00000000U,
                0x00000000U, 0xe9528d85U>;

using sfmt86243_64 =
    sfmt_engine<std::uint_fast64_t, 64, 86243, 366, 6, 56U, 19, 8U, 0xfdbffbffU,
                0xbff7ff3fU, 0xfd77efffU, 0xbf9ff3ffU, 0x00000001U, 0x00000000U,
                0x00000000U, 0xe9528d85U>;

using sfmt132049 =
    sfmt_engine<std::uint_fast32_t, 32, 132049, 110, 19, 8U, 21, 8U,
                0xffffbb5fU, 0xfb6ebf95U, 0xfffefffaU, 0xcff77fffU, 0x00000001U,
                0x00000000U, 0xcb520000U, 0xc7e91c7dU>;

using sfmt132049_64 =
    sfmt_engine<std::uint_fast64_t, 64, 132049, 110, 19, 8U, 21, 8U,
                0xffffbb5fU, 0xfb6ebf95U, 0xfffefffaU, 0xcff77fffU, 0x00000001U,
                0x00000000U, 0xcb520000U, 0xc7e91c7dU>;

using sfmt216091 =
    sfmt_engine<std::uint_fast32_t, 32, 216091, 627, 11, 24U, 10, 8U,
                0xbff7bff7U, 0xbfffffffU, 0xbffffa7fU, 0xffddfbfbU, 0xf8000001U,
                0x89e80709U, 0x3bd2b64bU, 0x0c64b1e4U>;

using sfmt216091_64 =
    sfmt_engine<std::uint_fast64_t, 64, 216091, 627, 11, 24U, 10, 8U,
                0xbff7bff7U, 0xbfffffffU, 0xbffffa7fU, 0xffddfbfbU, 0xf8000001U,
                0x89e80709U, 0x3bd2b64bU, 0x0c64b1e4U>;

} // namespace adhoc

#endif // BRNG_SFMT_ENGINE
