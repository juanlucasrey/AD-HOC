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

#ifndef BRNG_TRIVIUM_ENGINE
#define BRNG_TRIVIUM_ENGINE

#include "tools/mask.hpp"
#include "tools/unshift.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType> class trivium_engine final {
  private:
    void seed(const std::uint8_t *seed_and_iv, int length) {
        union SeedVec {
            std::uint8_t as8[16];
            UIntType as64[2];
        };
        SeedVec s;
        int elen = length > 10 ? 10 : length;
        for (int i = 0; i < 6; i++) {
            s.as8[i] = 0;
        }

        for (int i = 0; i < elen; i++) {
            s.as8[i + 6] = seed_and_iv[i];
        }

        for (int i = elen; i < 10; i++) {
            s.as8[i + 6] = 0;
        }

        for (int i = 0; i < 2; i++) {
            this->a[i] = s.as64[1 - i];
        }

        length -= 10;
        seed_and_iv += 10;
        length = length > 0 ? length : 0;
        SeedVec iv;
        for (int i = 0; i < 16 - length; i++) {
            iv.as8[i] = 0;
        }

        for (int i = 16 - length; i < 16; i++) {
            iv.as8[i] = seed_and_iv[i - (16 - length)];
        }

        for (int i = 0; i < 2; i++) {
            this->b[i] = iv.as64[1 - i];
        }

        c[0] = 0;
        c[1] = UIntType(7) << (128 - 111);

        for (int i = 0; i < 18; i++) {
            this->operator()<true>();
        }
    }

  public:
    using result_type = UIntType;
    static constexpr std::size_t word_size = 64U;
    static_assert(word_size <= std::numeric_limits<UIntType>::digits);
    static constexpr UIntType default_seed = 5489U;

    trivium_engine() : trivium_engine(default_seed) {}

    explicit trivium_engine(UIntType s) {
        std::uint8_t vec[8];

        for (int i = 0; i < 8; i++) {
            vec[i] = std::uint8_t(s);
            s >>= 8;
        }
        seed(vec, 8);
    }

    // template <class SeedSeq> explicit trivium_engine(SeedSeq &seq) {
    //     std::array<UIntType, 16> init_array{};
    //     seq.generate(init_array.begin(), init_array.end());
    //     // this->init(init_array);
    // }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        auto shift_array64 = [](std::array<UIntType, 3> const &vec,
                                std::size_t bits) -> UIntType {
            return (vec[0] << bits) | (vec[1] >> (64 - bits));
        };

        auto shift_array64_up = [](std::array<UIntType, 3> const &vec,
                                   std::size_t bits) -> UIntType {
            return (vec[0] << bits);
        };

        auto shift_array64_dn = [](std::array<UIntType, 3> const &vec,
                                   std::size_t bits) -> UIntType {
            return (vec[1] >> (64 - bits));
        };

        if constexpr (FwdDirection) {
            a[2] = shift_array64(c, 2) ^ shift_array64(c, 47);
            b[2] = shift_array64(a, 2) ^ shift_array64(a, 29);
            c[2] = shift_array64(b, 5) ^ shift_array64(b, 20);

            UIntType const new_a =
                a[2] ^ shift_array64(a, 5) ^
                (shift_array64(c, 46) & shift_array64(c, 45));
            UIntType const new_b =
                b[2] ^ shift_array64(b, 14) ^
                (shift_array64(a, 28) & shift_array64(a, 27));
            UIntType const new_c =
                c[2] ^ shift_array64(c, 23) ^
                (shift_array64(b, 19) & shift_array64(b, 18));

            a[1] = a[0];
            a[0] = new_a;
            b[1] = b[0];
            b[0] = new_b;
            c[1] = c[0];
            c[0] = new_c;

            return a[2] ^ b[2] ^ c[2];
        } else {
            UIntType const result = a[2] ^ b[2] ^ c[2];

            c[0] = c[1];
            b[0] = b[1];
            a[0] = a[1];

            c[1] = a[2] ^ (c[0] << 2) ^ (c[0] << 47);
            c[1] = unshift_right_xor<62, 45>(c[1]);
            // upper 47 are good
            c[1] <<= 17;

            a[1] = b[2] ^ (a[0] << 2) ^ (a[0] << 29);
            a[1] = unshift_right_xor<62, 27>(a[1]);
            // upper 29 are good
            a[1] <<= 35;

            b[1] = c[2] ^ (b[0] << 5) ^ (b[0] << 20);
            b[1] = unshift_right_xor<59, 15>(b[1]);
            // upper 20 are good
            b[1] <<= 44;

            std::array<UIntType, 3> a_p{};
            std::array<UIntType, 3> b_p{};
            std::array<UIntType, 3> c_p{};

            // upper 29 are good
            a_p[0] = a[1];

            // upper 20 are good
            b_p[0] = b[1];

            // upper 47 are good
            c_p[0] = c[1];

            UIntType const new_a_p = a[0];
            UIntType const new_b_p = b[0];
            UIntType const new_c_p = c[0];

            // update every 24 (29-5)
            UIntType a2_p_1 = 0;
            // update every 45 (47-2)
            UIntType a22_p_1 = 0;
            UIntType a2_xor = 0;

            // update every 6 (20-14)
            UIntType b2_p_1 = 0;
            // update every 27 (29-2)
            UIntType b22_p_1 = 0;
            UIntType b2_xor = 0;

            // update every 24 (47-23)
            UIntType c2_p_1 = 0;
            // update every 15 (20-5)
            UIntType c22_p_1 = 0;
            UIntType c2_xor = 0;

            std::array<bool, 2> b_bits{};
            constexpr auto m44 = mask<UIntType>(1, 44);
            b_bits[0] = m44 & b_p[0];
            constexpr auto m45 = mask<UIntType>(1, 45);
            b_bits[1] = m45 & b_p[0];

            std::array<bool, 2> a_bits{};
            constexpr auto m26 = mask<UIntType>(1, 35);
            a_bits[0] = m26 & a_p[0];
            constexpr auto m27 = mask<UIntType>(1, 36);
            a_bits[1] = m27 & a_p[0];

            std::array<bool, 2> c_bits{};
            constexpr auto m17 = mask<UIntType>(1, 17);
            c_bits[0] = m17 & c_p[0];
            constexpr auto m18 = mask<UIntType>(1, 18);
            c_bits[1] = m18 & c_p[0];

            constexpr bool templated = false;

            if constexpr (templated) {
                auto decrypt = [&]<std::size_t I>() {
                    constexpr std::size_t id1 = (I + 1) % 2;
                    constexpr std::size_t b_p_idx = I < 44 ? 0 : 1;
                    constexpr std::size_t a_p_idx = I < 35 ? 0 : 1;
                    constexpr std::size_t c_p_idx = I < 17 ? 0 : 1;
                    constexpr auto m = mask<UIntType>(1, 63 - I);
                    constexpr bool mult24 = ((I % 24) == 0);
                    constexpr bool mult45 = ((I % 45) == 0);
                    constexpr bool mult6 = ((I % 6) == 0);
                    constexpr bool mult27 = ((I % 27) == 0);
                    constexpr bool mult15 = ((I % 15) == 0);
                    constexpr auto mc = I < 17 ? mask<UIntType>(1, 16 - I)
                                               : mask<UIntType>(1, 63 + 17 - I);

                    constexpr auto ma = I < 35 ? mask<UIntType>(1, 34 - I)
                                               : mask<UIntType>(1, 63 + 35 - I);

                    constexpr auto mb = I < 44 ? mask<UIntType>(1, 43 - I)
                                               : mask<UIntType>(1, 63 + 44 - I);

                    if constexpr (mult24) {
                        if constexpr (I + 24 < 59) {
                            a2_p_1 = new_a_p ^ shift_array64_up(a_p, 5);
                        } else if constexpr (I < 59) {
                            a2_p_1 = new_a_p ^ shift_array64(a_p, 5);
                        } else {
                            a2_p_1 = new_a_p ^ shift_array64_dn(a_p, 5);
                        }
                    }

                    if constexpr (mult45) {
                        if constexpr (I + 45 < 62) {
                            a22_p_1 = shift_array64_up(c_p, 2);
                        } else if constexpr (I < 62) {
                            a22_p_1 = shift_array64(c_p, 2);
                        } else {
                            a22_p_1 = shift_array64_dn(c_p, 2);
                        }
                    }

                    if constexpr (mult24 || mult45) {
                        a2_xor = a2_p_1 ^ a22_p_1;
                    }

                    if constexpr (mult6) {
                        if constexpr (I + 6 < 50) {
                            b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                        } else if constexpr (I < 50) {
                            b2_p_1 = new_b_p ^ shift_array64(b_p, 14);
                        } else {
                            b2_p_1 = new_b_p ^ shift_array64_dn(b_p, 14);
                        }
                    }

                    if constexpr (mult27) {
                        if constexpr (I + 27 < 62) {
                            b22_p_1 = shift_array64_up(a_p, 2);
                        } else if constexpr (I < 62) {
                            b22_p_1 = shift_array64(a_p, 2);
                        } else {
                            b22_p_1 = shift_array64_dn(a_p, 2);
                        }
                    }

                    if constexpr (mult6 || mult27) {
                        b2_xor = b2_p_1 ^ b22_p_1;
                    }

                    if constexpr (mult24) {
                        if constexpr (I + 24 < 41) {
                            c2_p_1 = new_c_p ^ shift_array64_up(c_p, 23);
                        } else if constexpr (I < 41) {
                            c2_p_1 = new_c_p ^ shift_array64(c_p, 23);
                        } else {
                            c2_p_1 = new_c_p ^ shift_array64_dn(c_p, 23);
                        }
                    }

                    if constexpr (mult15) {
                        if constexpr (I + 15 < 59) {
                            c22_p_1 = shift_array64_up(b_p, 5);
                        } else if constexpr (I < 59) {
                            c22_p_1 = shift_array64(b_p, 5);
                        } else {
                            c22_p_1 = shift_array64_dn(b_p, 5);
                        }
                    }

                    if constexpr (mult24 || mult15) {
                        c2_xor = c2_p_1 ^ c22_p_1;
                    }

                    b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                  (b_bits[0] && b_bits[1]);
                    b_p[b_p_idx] |= (mb * b_bits[id1]);

                    a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                  (a_bits[0] && a_bits[1]);
                    a_p[a_p_idx] |= (ma * a_bits[id1]);

                    c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                  (c_bits[0] && c_bits[1]);
                    c_p[c_p_idx] |= (mc * c_bits[id1]);
                };

                // cannot call a lambda template recursively..
                decrypt.template operator()<0>();
                decrypt.template operator()<1>();
                decrypt.template operator()<2>();
                decrypt.template operator()<3>();
                decrypt.template operator()<4>();
                decrypt.template operator()<5>();
                decrypt.template operator()<6>();
                decrypt.template operator()<7>();
                decrypt.template operator()<8>();
                decrypt.template operator()<9>();
                decrypt.template operator()<10>();
                decrypt.template operator()<11>();
                decrypt.template operator()<12>();
                decrypt.template operator()<13>();
                decrypt.template operator()<14>();
                decrypt.template operator()<15>();
                decrypt.template operator()<16>();
                decrypt.template operator()<17>();
                decrypt.template operator()<18>();
                decrypt.template operator()<19>();
                decrypt.template operator()<20>();
                decrypt.template operator()<21>();
                decrypt.template operator()<22>();
                decrypt.template operator()<23>();
                decrypt.template operator()<24>();
                decrypt.template operator()<25>();
                decrypt.template operator()<26>();
                decrypt.template operator()<27>();
                decrypt.template operator()<28>();
                decrypt.template operator()<29>();
                decrypt.template operator()<30>();
                decrypt.template operator()<31>();
                decrypt.template operator()<32>();
                decrypt.template operator()<33>();
                decrypt.template operator()<34>();
                decrypt.template operator()<35>();
                decrypt.template operator()<36>();
                decrypt.template operator()<37>();
                decrypt.template operator()<38>();
                decrypt.template operator()<39>();
                decrypt.template operator()<40>();
                decrypt.template operator()<41>();
                decrypt.template operator()<42>();
                decrypt.template operator()<43>();
                decrypt.template operator()<44>();
                decrypt.template operator()<45>();
                decrypt.template operator()<46>();
                decrypt.template operator()<47>();
                decrypt.template operator()<48>();
                decrypt.template operator()<49>();
                decrypt.template operator()<50>();
                decrypt.template operator()<51>();
                decrypt.template operator()<52>();
                decrypt.template operator()<53>();
                decrypt.template operator()<54>();
                decrypt.template operator()<55>();
                decrypt.template operator()<56>();
                decrypt.template operator()<57>();
                decrypt.template operator()<58>();
                decrypt.template operator()<59>();
                decrypt.template operator()<60>();
                decrypt.template operator()<61>();
                decrypt.template operator()<62>();
                decrypt.template operator()<63>();

            } else {

                a2_p_1 = new_a_p ^ shift_array64_up(a_p, 5);
                a22_p_1 = shift_array64_up(c_p, 2);
                a2_xor = a2_p_1 ^ a22_p_1;
                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b22_p_1 = shift_array64_up(a_p, 2);
                b2_xor = b2_p_1 ^ b22_p_1;
                c2_p_1 = new_c_p ^ shift_array64_up(c_p, 23);
                c22_p_1 = shift_array64_up(b_p, 5);
                c2_xor = c2_p_1 ^ c22_p_1;

                auto m = mask<UIntType>(1, 63);
                auto mb = mask<UIntType>(1, 43);
                auto ma = mask<UIntType>(1, 34);
                auto mc = mask<UIntType>(1, 16);

                std::size_t i = 0;
                for (; i < 6; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[0] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;

                for (; i < 12; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);
                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[0] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;

                for (; i < 15; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[0] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                c22_p_1 = shift_array64_up(b_p, 5);
                c2_xor = c2_p_1 ^ c22_p_1;

                for (; i < 17; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);
                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[0] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                mc = mask<UIntType>(1, 63);

                for (; i < 18; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);
                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;

                for (; i < 19; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);
                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                for (; i < 24; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                a2_p_1 = new_a_p ^ shift_array64_up(a_p, 5);
                a2_xor = a2_p_1 ^ a22_p_1;
                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;
                c2_p_1 = new_c_p ^ shift_array64(c_p, 23);
                c2_xor = c2_p_1 ^ c22_p_1;

                for (; i < 27; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);
                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b22_p_1 = shift_array64_up(a_p, 2);
                b2_xor = b2_p_1 ^ b22_p_1;

                for (; i < 30; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;
                c22_p_1 = shift_array64_up(b_p, 5);
                c2_xor = c2_p_1 ^ c22_p_1;

                for (; i < 35; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[0] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                ma = mask<UIntType>(1, 63);

                for (; i < 36; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;

                for (; i < 37; ++i) {
                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                for (; i < 42; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b2_p_1 = new_b_p ^ shift_array64_up(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;

                for (; i < 44; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);
                        b_p[0] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                mb = mask<UIntType>(1, 63);

                for (; i < 45; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);
                        b_p[1] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                a22_p_1 = shift_array64(c_p, 2);
                a2_xor = a2_p_1 ^ a22_p_1;
                c22_p_1 = shift_array64(b_p, 5);
                c2_xor = c2_p_1 ^ c22_p_1;

                for (; i < 46; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[1] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                for (; i < 48; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[1] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                a2_p_1 = new_a_p ^ shift_array64(a_p, 5);
                a2_xor = a2_p_1 ^ a22_p_1;
                b2_p_1 = new_b_p ^ shift_array64(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;
                c2_p_1 = new_c_p ^ shift_array64_dn(c_p, 23);
                c2_xor = c2_p_1 ^ c22_p_1;

                for (; i < 54; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[1] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b22_p_1 = shift_array64(a_p, 2);
                b2_p_1 = new_b_p ^ shift_array64_dn(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;

                for (; i < 60; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[1] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }

                b2_p_1 = new_b_p ^ shift_array64_dn(b_p, 14);
                b2_xor = b2_p_1 ^ b22_p_1;
                c22_p_1 = shift_array64_dn(b_p, 5);
                c2_xor = c2_p_1 ^ c22_p_1;

                for (; i < 64; ++i) {

                    std::size_t id1 = (i + 1) % 2;
                    {
                        b_bits[id1] = static_cast<bool>(c2_xor & m) !=
                                      (b_bits[0] && b_bits[1]);

                        b_p[1] |= (mb * b_bits[id1]);

                        mb >>= 1;
                    }

                    {
                        a_bits[id1] = static_cast<bool>(b2_xor & m) !=
                                      (a_bits[0] && a_bits[1]);

                        a_p[1] |= (ma * a_bits[id1]);

                        ma >>= 1;
                    }

                    {
                        c_bits[id1] = static_cast<bool>(a2_xor & m) !=
                                      (c_bits[0] && c_bits[1]);

                        c_p[1] |= (mc * c_bits[id1]);

                        mc >>= 1;
                    }

                    m >>= 1;
                }
            }

            a[1] = a_p[0];
            b[1] = b_p[0];
            c[1] = c_p[0];

            a[2] = shift_array64(c_p, 2) ^ shift_array64(c_p, 47);
            b[2] = shift_array64(a_p, 2) ^ shift_array64(a_p, 29);
            c[2] = shift_array64(b_p, 5) ^ shift_array64(b_p, 20);

            return result;
        }
    }

    template <bool FwdDirection = true> void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()<FwdDirection>();
        }
    }

    static constexpr auto min() -> result_type {
        return static_cast<result_type>(0U);
    }

    static constexpr auto max() -> result_type {
        return mask<UIntType>(word_size);
    }

    auto operator==(const trivium_engine &rhs) const -> bool {
        return (this->a == rhs.a) && (this->b == rhs.b) && (this->c == rhs.c);
    }

    auto operator!=(const trivium_engine &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    std::array<UIntType, 3> a{};
    std::array<UIntType, 3> b{};
    std::array<UIntType, 3> c{};
};

} // namespace adhoc

#endif // BRNG_TRIVIUM_ENGINE
