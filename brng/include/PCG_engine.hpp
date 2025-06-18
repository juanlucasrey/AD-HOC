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

#ifndef BRNG_PCG_ENGINE
#define BRNG_PCG_ENGINE

#include "tools/modular_arithmetic.hpp"
#include "tools/uint128.hpp"

#include <array>
#include <cstdint>

namespace adhoc {

namespace detail {

template <class UIntType, std::size_t w, class UpgradedType>
auto constexpr init_upgraded(UIntType low, UIntType high) -> UpgradedType {
    UpgradedType const low_up = static_cast<UpgradedType>(low);
    UpgradedType const high_up = static_cast<UpgradedType>(high);
    return low_up + (high_up << w);
}

using bitcount_t = std::uint16_t;

template <typename itype>
inline auto rotr(itype value, bitcount_t rot) -> itype {
    constexpr bitcount_t bits = sizeof(itype) * 8;
    constexpr bitcount_t mask = bits - 1;
    return (value >> rot) | (value << ((-rot) & mask));
}

template <typename xtype, typename itype> auto xsh_rr(itype internal) -> xtype {
    constexpr auto bits = bitcount_t{sizeof(itype) * 8};
    constexpr auto xtypebits = bitcount_t{sizeof(xtype) * 8};
    constexpr bitcount_t sparebits = bits - xtypebits;
    constexpr bitcount_t wantedopbits = xtypebits >= 128  ? 7
                                        : xtypebits >= 64 ? 6
                                        : xtypebits >= 32 ? 5
                                        : xtypebits >= 16 ? 4
                                                          : 3;
    constexpr bitcount_t opbits =
        sparebits >= wantedopbits ? wantedopbits : sparebits;
    constexpr bitcount_t amplifier = wantedopbits - opbits;
    constexpr bitcount_t mask = (1U << opbits) - 1U;
    constexpr bitcount_t topspare = opbits;
    constexpr bitcount_t bottomspare = sparebits - topspare;
    constexpr bitcount_t xshift = (topspare + xtypebits) / 2;

    bitcount_t rot =
        opbits ? bitcount_t(internal >> (bits - opbits)) & mask : 0;
    bitcount_t amprot = (rot << amplifier) & mask;
    internal ^= internal >> xshift;
    xtype result = xtype(internal >> bottomspare);
    result = rotr(result, amprot);
    return result;
}

template <typename xtype, typename itype> auto xsh_rs(itype internal) -> xtype {
    constexpr auto bits = bitcount_t(sizeof(itype) * 8);
    constexpr auto xtypebits = bitcount_t(sizeof(xtype) * 8);
    constexpr bitcount_t sparebits = bits - xtypebits;
    constexpr bitcount_t opbits = sparebits - 5 >= 64   ? 5
                                  : sparebits - 4 >= 32 ? 4
                                  : sparebits - 3 >= 16 ? 3
                                  : sparebits - 2 >= 4  ? 2
                                  : sparebits - 1 >= 1  ? 1
                                                        : 0;
    constexpr bitcount_t mask = (1U << opbits) - 1;
    constexpr bitcount_t maxrandshift = mask;
    constexpr bitcount_t topspare = opbits;
    constexpr bitcount_t bottomspare = sparebits - topspare;
    constexpr bitcount_t xshift = topspare + ((xtypebits + maxrandshift) / 2);
    bitcount_t rshift =
        opbits ? bitcount_t(internal >> (bits - opbits)) & mask : 0;
    internal ^= internal >> xshift;
    xtype result = xtype(internal >> (bottomspare - maxrandshift + rshift));
    return result;
}

template <typename xtype, typename itype> auto xsl_rr(itype internal) -> xtype {
    constexpr bitcount_t xtypebits = bitcount_t(sizeof(xtype) * 8);
    constexpr bitcount_t bits = bitcount_t(sizeof(itype) * 8);
    constexpr bitcount_t sparebits = bits - xtypebits;
    constexpr bitcount_t wantedopbits = xtypebits >= 128  ? 7
                                        : xtypebits >= 64 ? 6
                                        : xtypebits >= 32 ? 5
                                        : xtypebits >= 16 ? 4
                                                          : 3;
    constexpr bitcount_t opbits =
        sparebits >= wantedopbits ? wantedopbits : sparebits;
    constexpr bitcount_t amplifier = wantedopbits - opbits;
    constexpr bitcount_t mask = (1 << opbits) - 1;
    constexpr bitcount_t topspare = sparebits;
    constexpr bitcount_t bottomspare = sparebits - topspare;
    constexpr bitcount_t xshift = (topspare + xtypebits) / 2;

    bitcount_t rot =
        opbits ? bitcount_t(internal >> (bits - opbits)) & mask : 0;
    bitcount_t amprot = (rot << amplifier) & mask;
    internal ^= internal >> xshift;
    xtype result = xtype(internal >> bottomspare);
    result = rotr(result, amprot);
    return result;
}

template <typename xtype, std::size_t w, typename itype>
auto rxs_m_xs(itype internal) -> xtype {
    constexpr bitcount_t xtypebits = bitcount_t(sizeof(xtype) * 8);
    constexpr bitcount_t bits = bitcount_t(sizeof(itype) * 8);
    constexpr bitcount_t opbits = xtypebits >= 128  ? 6
                                  : xtypebits >= 64 ? 5
                                  : xtypebits >= 32 ? 4
                                  : xtypebits >= 16 ? 3
                                                    : 2;
    constexpr bitcount_t shift = bits - xtypebits;
    constexpr bitcount_t mask = (1 << opbits) - 1;
    bitcount_t rshift =
        opbits ? bitcount_t(internal >> (bits - opbits)) & mask : 0;
    internal ^= internal >> (opbits + rshift);

    if constexpr (w <= 32) {
        internal *= 277803737U;
    } else if constexpr (w <= 64) {
        internal *= 12605985483714917081ULL;
    } else if constexpr (w <= 128) {
        constexpr itype multiplier = init_upgraded<std::uint64_t, 64, uint128>(
            12605985483714917081ULL, 17766728186571221404ULL);
        internal *= multiplier;
    }

    xtype result = internal >> shift;
    result ^= result >> ((2U * xtypebits + 2U) / 3U);
    return result;
}

template <typename itype>
inline auto unxorshift(itype x, bitcount_t bits, bitcount_t shift) -> itype {
    if (2 * shift >= bits) {
        return x ^ (x >> shift);
    }
    itype lowmask1 = (itype(1U) << (bits - shift * 2)) - 1;
    itype highmask1 = ~lowmask1;
    itype top1 = x;
    itype bottom1 = x & lowmask1;
    top1 ^= top1 >> shift;
    top1 &= highmask1;
    x = top1 | bottom1;
    itype lowmask2 = (itype(1U) << (bits - shift)) - 1;
    itype bottom2 = x & lowmask2;
    bottom2 = unxorshift(bottom2, bits - shift, shift);
    bottom2 &= lowmask1;
    return top1 | bottom2;
}

template <std::size_t w, typename itype>
auto rxs_m_xs_back(itype internal) -> itype {
    constexpr bitcount_t bits = bitcount_t(sizeof(itype) * 8);
    constexpr bitcount_t opbits = bits >= 128  ? 6
                                  : bits >= 64 ? 5
                                  : bits >= 32 ? 4
                                  : bits >= 16 ? 3
                                               : 2;
    constexpr bitcount_t mask = (1 << opbits) - 1;

    internal = unxorshift(internal, bits, (2U * bits + 2U) / 3U);

    // internal *= mcg_unmultiplier<itype>::unmultiplier();
    if constexpr (w <= 32) {
        internal *= 2897767785U;
    } else if constexpr (w <= 64) {
        internal *= 15009553638781119849ULL;
    } else if constexpr (w <= 128) {
        constexpr itype multiplier = init_upgraded<std::uint64_t, 64, uint128>(
            15009553638781119849ULL, 14422606686972528997ULL);
        internal *= multiplier;
    }

    bitcount_t rshift = opbits ? (internal >> (bits - opbits)) & mask : 0;
    internal = unxorshift(internal, bits, opbits + rshift);

    return internal;
}

template <typename T> struct halfsize_trait {};
template <> struct halfsize_trait<uint128> {
    using type = uint64_t;
};
template <> struct halfsize_trait<uint64_t> {
    using type = uint32_t;
};
template <> struct halfsize_trait<uint32_t> {
    using type = uint16_t;
};
template <> struct halfsize_trait<uint16_t> {
    using type = uint8_t;
};

template <typename xtype, typename itype>
auto xsl_rr_rr(itype internal) -> itype {
    using htype = typename halfsize_trait<itype>::type;

    constexpr bitcount_t htypebits = bitcount_t(sizeof(htype) * 8);
    constexpr bitcount_t bits = bitcount_t(sizeof(itype) * 8);
    constexpr bitcount_t sparebits = bits - htypebits;
    constexpr bitcount_t wantedopbits = htypebits >= 128  ? 7
                                        : htypebits >= 64 ? 6
                                        : htypebits >= 32 ? 5
                                        : htypebits >= 16 ? 4
                                                          : 3;
    constexpr bitcount_t opbits =
        sparebits >= wantedopbits ? wantedopbits : sparebits;
    constexpr bitcount_t amplifier = wantedopbits - opbits;
    constexpr bitcount_t mask = (1 << opbits) - 1;
    constexpr bitcount_t topspare = sparebits;
    constexpr bitcount_t xshift = (topspare + htypebits) / 2;

    bitcount_t rot =
        opbits ? bitcount_t(internal >> (bits - opbits)) & mask : 0;
    bitcount_t amprot = (rot << amplifier) & mask;
    internal ^= internal >> xshift;
    htype lowbits = htype(internal);
    lowbits = rotr(lowbits, amprot);
    htype highbits = htype(internal >> topspare);
    bitcount_t rot2 = lowbits & mask;
    bitcount_t amprot2 = (rot2 << amplifier) & mask;
    highbits = rotr(highbits, amprot2);
    return (itype(highbits) << topspare) ^ itype(lowbits);
}

} // namespace detail

enum class tempering_type {
    xsh_rr,
    xsh_rs,
    xsl_rr,
    rxs_m_xs,
    xsl_rr_rr,
};

template <class UIntType = std::uint32_t, std::size_t w = 32U,
          tempering_type TemperingType = tempering_type::xsh_rr,
          bool upgrade = true, std::size_t table_pow2 = 0U,
          std::size_t advance_pow2 = 0U, bool upgrade_table = false,
          bool mcg = false, bool kdd = true, bool invert_upgrade_call = false>
class PCG_engine final {

    static_assert(w <= std::numeric_limits<UIntType>::digits);

    static constexpr std::size_t w_upgraded =
        upgrade ? ((w <= 32U) ? 64U : 128U) : w;

    static_assert(w_upgraded <= 128U);

    static constexpr std::size_t w_table_upgraded =
        upgrade_table ? ((w <= 32U) ? 64U : 128U) : w;

    static_assert(w_table_upgraded <= 128U);

    using upgraded_type = std::conditional_t<
        w_upgraded == 32, std::uint32_t,
        std::conditional_t<w_upgraded == 64, std::uint64_t, uint128>>;

    using upgraded_type_table = std::conditional_t<
        w_table_upgraded == 32, std::uint32_t,
        std::conditional_t<w_table_upgraded == 64, std::uint64_t, uint128>>;

    static constexpr upgraded_type multiplier =
        (w_upgraded == 32) ? upgraded_type(747796405U)
        : (w_upgraded == 64)
            ? upgraded_type(6364136223846793005ULL)
            : upgraded_type(detail::init_upgraded<std::uint64_t, 64, uint128>(
                  4865540595714422341ULL, 2549297995355413924ULL));

    static constexpr upgraded_type_table multiplier_table =
        (w_table_upgraded == 32) ? upgraded_type_table(747796405U)
        : (w_table_upgraded == 64)
            ? upgraded_type_table(6364136223846793005ULL)
            : upgraded_type_table(
                  detail::init_upgraded<std::uint64_t, 64, uint128>(
                      4865540595714422341ULL, 2549297995355413924ULL));

    static constexpr upgraded_type_table increment_table =
        (w_table_upgraded == 32) ? upgraded_type_table(2891336453U)
        : (w_table_upgraded == 64)
            ? upgraded_type_table(1442695040888963407ULL)
            : upgraded_type_table(
                  detail::init_upgraded<std::uint64_t, 64, uint128>(
                      1442695040888963407ULL, 6364136223846793005ULL));

    static constexpr std::size_t table_size = 1UL << table_pow2;
    static constexpr upgraded_type table_mask =
        (upgraded_type(1U) << table_pow2) - upgraded_type(1U);

    static constexpr std::size_t stypebits = sizeof(upgraded_type) * 8;
    static constexpr std::size_t table_shift = stypebits - table_pow2;

    static constexpr std::size_t tick_limit_pow2 = 64U;

    void selfinit() {
        if constexpr (table_pow2) {
            result_type lhs = this->sub_operator();
            result_type rhs = this->sub_operator();
            result_type xdiff = lhs - rhs;
            for (std::size_t i = 0; i < table_size; ++i) {
                data_[i] = this->sub_operator() ^ xdiff;
            }
        }
    }

  public:
    using result_type = UIntType;

    static constexpr UIntType default_state_low =
        w <= 32 ? 3505598835U : 10954001071364377171ULL;
    static constexpr UIntType default_state_high =
        w <= 32 ? 1297702388U : 13320540357896250258ULL;

    static constexpr UIntType default_seq_low =
        w <= 32 ? 4150755663U : 1442695040888963407ULL;
    static constexpr UIntType default_seq_high =
        w <= 32 ? 335903614U : 6364136223846793005ULL;

    PCG_engine()
        : PCG_engine(default_state_low, default_state_high, default_seq_low,
                     default_seq_high) {}

    explicit PCG_engine(upgraded_type initstate, upgraded_type initseq)
        : state(initstate), inc(initseq) {
        this->selfinit();
    }

    explicit PCG_engine(UIntType initstate_low, UIntType initstate_high,
                        UIntType initseq_low, UIntType initseq_high)
        : state(detail::init_upgraded<UIntType, w, upgraded_type>(
              initstate_low, initstate_high)),
          inc(detail::init_upgraded<UIntType, w, upgraded_type>(initseq_low,
                                                                initseq_high)) {
        static_assert(upgrade);
        this->selfinit();
    }

    template <bool FwdDirection = true>
    auto get_extended_value() -> result_type {
        upgraded_type state_temp = this->state;
        if (kdd && mcg) {
            // The low order bits of an MCG are constant, so drop them.
            state_temp >>= 2;
        }

        std::size_t index =
            kdd ? state_temp & table_mask : state_temp >> table_shift;

        result_type result = 0;
        if constexpr (!FwdDirection) {
            result = data_[index];
        }

        static constexpr bool may_tick =
            (advance_pow2 < stypebits) && (advance_pow2 < tick_limit_pow2);
        if constexpr (may_tick) {
            static constexpr std::size_t tick_shift = stypebits - advance_pow2;
            static constexpr upgraded_type tick_mask =
                upgraded_type((std::uint64_t(1) << advance_pow2) - 1);

            bool tick = kdd ? (state_temp & tick_mask) == upgraded_type(0u)
                            : (state_temp >> tick_shift) == upgraded_type(0u);
            if (tick) {
                bool carry = false;
                for (std::size_t i = 0; i < table_size; ++i) {
                    upgraded_type_table new_incr =
                        upgraded_type_table((i + 1) * 2);

                    bool new_carry;

                    upgraded_type_table state =
                        detail::rxs_m_xs_back<w>(this->data_[i]);

                    if constexpr (FwdDirection) {
                        state = state * multiplier_table + increment_table +
                                new_incr;
                        this->data_[i] =
                            detail::rxs_m_xs<result_type, w>(state);
                        upgraded_type_table const zero =
                            mcg ? state & upgraded_type_table(3U)
                                : upgraded_type_table(0U);

                        new_carry = this->data_[i] == zero;
                    } else {
                        upgraded_type_table const zero =
                            mcg ? state & upgraded_type_table(3U)
                                : upgraded_type_table(0U);

                        new_carry = this->data_[i] == zero;

                        state -= increment_table + new_incr;

                        static constexpr upgraded_type
                            multiplier_table_inverse =
                                modular_multiplicative_inverse_pow2(
                                    w_table_upgraded, multiplier_table);
                        state *= multiplier_table_inverse;

                        this->data_[i] =
                            detail::rxs_m_xs<result_type, w>(state);
                    }

                    if (carry) {
                        carry = new_carry;
                    } else {
                        carry = carry || new_carry;
                    }
                }
            }
        }

        if constexpr (FwdDirection) {
            return data_[index];
        } else {
            return result;
        }
    }

    template <bool FwdDirection = true>
    inline auto sub_operator() -> result_type {
        upgraded_type result = 0;
        if constexpr (FwdDirection) {
            this->state = this->state * multiplier + this->inc;
            result = this->state;
        } else {
            static constexpr upgraded_type multiplier_inverse =
                modular_multiplicative_inverse_pow2(w_upgraded, multiplier);

            result = this->state;
            this->state = (this->state - this->inc) * multiplier_inverse;
        }

        if constexpr (TemperingType == tempering_type::xsh_rr) {
            return detail::xsh_rr<result_type>(result);
        } else if constexpr (TemperingType == tempering_type::xsh_rs) {
            return detail::xsh_rs<result_type>(result);
        } else if constexpr (TemperingType == tempering_type::xsl_rr) {
            return detail::xsl_rr<result_type>(result);
        } else if constexpr (TemperingType == tempering_type::rxs_m_xs) {
            return detail::rxs_m_xs<result_type, w>(result);
        } else if constexpr (TemperingType == tempering_type::xsl_rr_rr) {
            return detail::xsl_rr_rr<result_type>(result);
        }
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (table_pow2) {
            if constexpr (FwdDirection == invert_upgrade_call) {
                result_type const lhs = this->sub_operator<FwdDirection>();
                result_type const rhs =
                    this->get_extended_value<FwdDirection>();
                return lhs ^ rhs;
            } else {
                result_type const rhs =
                    this->get_extended_value<FwdDirection>();
                result_type const lhs = this->sub_operator<FwdDirection>();
                return lhs ^ rhs;
            }
        } else {
            return this->sub_operator<FwdDirection>();
        }
    }

    auto operator==(const PCG_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->inc == rhs.inc) &&
               (this->data_ == rhs.data_);
        ;
    }

    auto operator!=(const PCG_engine &rhs) -> bool {
        return !(this->operator==(rhs));
    }

  private:
    upgraded_type state;
    upgraded_type inc;

    std::array<UIntType, table_size> data_{};
};

} // namespace adhoc

#endif // BRNG_PCG_ENGINE
