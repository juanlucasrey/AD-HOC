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

#ifndef BRNG_SOBOL_ENGINE
#define BRNG_SOBOL_ENGINE

#include "sobol_tables/joe_kuo_old_1111.hpp"
#include "sobol_tables/joe_kuo_other_0_7600.hpp"
#include "sobol_tables/joe_kuo_other_2_3900.hpp"
#include "sobol_tables/joe_kuo_other_3_7300.hpp"
#include "sobol_tables/joe_kuo_other_4_5600.hpp"
#include "sobol_tables/new_joe_kuo_5_21201.hpp"
#include "sobol_tables/new_joe_kuo_6_21201.hpp"
#include "sobol_tables/new_joe_kuo_7_21201.hpp"
#include "tools/bit_find.hpp"
#include "tools/mask.hpp"
#include "tools/uint128.hpp"

#include <array>
#include <bitset>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace adhoc {

enum class qrng_table {
    joe_kuo_old_1111, // a.k.a. joe_kuo_2003.
    joe_kuo_other_0_7600,
    joe_kuo_other_2_3900,
    joe_kuo_other_3_7300,
    joe_kuo_other_4_5600,
    new_joe_kuo_5_21201,
    new_joe_kuo_6_21201, // a.k.a. joe_kuo_2008. recommended.
    new_joe_kuo_7_21201
};

namespace detail {

constexpr auto max_dimension(qrng_table table) -> unsigned int {
    switch (table) {
    case qrng_table::joe_kuo_old_1111:
        return joe_kuo_old_1111_poly.size();
    case qrng_table::joe_kuo_other_0_7600:
        return joe_kuo_other_0_7600_poly.size();
    case qrng_table::joe_kuo_other_2_3900:
        return joe_kuo_other_2_3900_poly.size();
    case qrng_table::joe_kuo_other_3_7300:
        return joe_kuo_other_3_7300_poly.size();
    case qrng_table::joe_kuo_other_4_5600:
        return joe_kuo_other_4_5600_poly.size();
    case qrng_table::new_joe_kuo_5_21201:
        return new_joe_kuo_5_21201_poly.size();
    case qrng_table::new_joe_kuo_6_21201:
        return new_joe_kuo_6_21201_poly.size();
    case qrng_table::new_joe_kuo_7_21201:
        return new_joe_kuo_7_21201_poly.size();
    }
}

inline auto poly(qrng_table table) {
    switch (table) {
    case qrng_table::joe_kuo_old_1111:
        return joe_kuo_old_1111_poly.begin();
    case qrng_table::joe_kuo_other_0_7600:
        return joe_kuo_other_0_7600_poly.begin();
    case qrng_table::joe_kuo_other_2_3900:
        return joe_kuo_other_2_3900_poly.begin();
    case qrng_table::joe_kuo_other_3_7300:
        return joe_kuo_other_3_7300_poly.begin();
    case qrng_table::joe_kuo_other_4_5600:
        return joe_kuo_other_4_5600_poly.begin();
    case qrng_table::new_joe_kuo_5_21201:
        return new_joe_kuo_5_21201_poly.begin();
    case qrng_table::new_joe_kuo_6_21201:
        return new_joe_kuo_6_21201_poly.begin();
    case qrng_table::new_joe_kuo_7_21201:
        return new_joe_kuo_7_21201_poly.begin();
    }
    return joe_kuo_old_1111_poly.begin();
}

inline auto init(qrng_table table) {
    switch (table) {
    case qrng_table::joe_kuo_old_1111:
        return joe_kuo_old_1111_init.begin();
    case qrng_table::joe_kuo_other_0_7600:
        return joe_kuo_other_0_7600_init.begin();
    case qrng_table::joe_kuo_other_2_3900:
        return joe_kuo_other_2_3900_init.begin();
    case qrng_table::joe_kuo_other_3_7300:
        return joe_kuo_other_3_7300_init.begin();
    case qrng_table::joe_kuo_other_4_5600:
        return joe_kuo_other_4_5600_init.begin();
    case qrng_table::new_joe_kuo_5_21201:
        return new_joe_kuo_5_21201_init.begin();
    case qrng_table::new_joe_kuo_6_21201:
        return new_joe_kuo_6_21201_init.begin();
    case qrng_table::new_joe_kuo_7_21201:
        return new_joe_kuo_7_21201_init.begin();
    }
    return joe_kuo_old_1111_init.begin();
}

inline auto mask(std::size_t w) -> std::uint32_t {
    if (w == 0) {
        return 0;
    }
    std::uint32_t ach = 1U << (w - 1U);

    return ach | (ach - 1U);

    // equivalent to the following but shift does not overflow
    // return (static_cast<result_type>(1U) << w) - 1U;
}

} // namespace detail

template <class UIntType, std::size_t w, bool skip_first = true,
          qrng_table t = qrng_table::joe_kuo_old_1111>
class sobol_engine final {
    static_assert(w <= std::numeric_limits<UIntType>::digits);
    static_assert(w <= 64);

    inline void increase_counter() {
        constexpr auto mask = sobol_engine::max();
        this->m_seq_num = (this->m_seq_num + 1) & mask;
    }

    inline void decrease_counter() {
        constexpr auto mask = sobol_engine::max();
        this->m_seq_num = this->m_seq_num == 0 ? mask : (this->m_seq_num - 1);
    }

    inline void generate() {
        constexpr auto digits = std::numeric_limits<UIntType>::digits;
        using ceiling_type =
            std::conditional_t<digits <= 32, std::uint32_t, std::uint64_t>;

        if constexpr (skip_first) {
            auto const c = this->m_seq_num
                               ? findLeastSignificantZero(
                                     static_cast<ceiling_type>(this->m_seq_num))
                               : (w + 1);

            for (std::size_t i = 0; i < this->N; ++i) {
                this->Y[i] ^= this->v[i][c];
            }
        } else {
            auto const c = findLeastSignificantZero(
                static_cast<ceiling_type>(this->m_seq_num));

            for (std::size_t i = 0; i < this->N; ++i) {
                this->Y[i] ^= this->v[i][c];
            }
        }
    }

  public:
    using result_type = UIntType;

    explicit sobol_engine(unsigned int n_dims = 1U)
        : N{n_dims}, v(n_dims), Y(n_dims) {

        constexpr unsigned int max_dim = detail::max_dimension(t);
        if (n_dims < 1 || n_dims > max_dim) {
            throw std::runtime_error("dimension error");
        }

        this->v.front().fill(1U);

        const auto *poly_p = detail::poly(t);
        poly_p++;

        const auto *init_p = detail::init(t);
        std::size_t idx = 0;
        std::uint32_t val = (*init_p++);

        for (std::size_t i = 1; i < this->N; ++i) {
            unsigned int const poly_e = (*poly_p++);
            std::size_t const poly_d = highestBitIndex(poly_e);
            std::bitset<w> includ(poly_e);

            std::size_t ii = 0;
            for (; ii < poly_d; ++ii) {
                std::uint32_t init = 0;
                {
                    std::size_t bits = ii;
                    std::size_t shift = 0;
                    std::size_t rembits = 32U - idx;
                    while (bits) {
                        std::size_t bits_here = std::min(rembits, bits);
                        std::uint32_t m = detail::mask(bits_here);
                        init |= (val & m) << shift;
                        val >>= bits_here;
                        idx += bits_here;
                        bits -= bits_here;
                        shift += bits_here;

                        if (idx == 32) {
                            idx = 0;
                            val = (*init_p++);
                        }
                        rembits = 32U - idx;
                    }

                    init = (init << 1U) | 1U;
                }

                this->v[i][ii] = init;
            }

            for (; ii < w; ++ii) {
                auto &v_new = this->v[i][ii];
                v_new = this->v[i][ii - poly_d];

                UIntType kk = 1;
                for (std::size_t k = 0; k < poly_d; ++k) {
                    kk *= 2;
                    if (includ[poly_d - k - 1]) {
                        v_new ^= (kk * v[i][ii - k - 1]);
                    }
                }
            }
        }

        UIntType kk = 1;
        for (std::size_t i = w - 1; i-- > 0;) {
            kk *= 2;

            for (std::size_t ii = 0; ii < this->N; ++ii) {
                this->v[ii][i] *= kk;
            }
        }

        if constexpr (skip_first) {
            for (std::size_t i = 0; i < this->N; ++i) {
                this->Y[i] = this->v[i][0];
            }

            for (std::size_t i = 0; i < this->N; ++i) {
                this->v[i][w] = this->v[i][w - 1] ^ this->v[i][0];
            }

            for (std::size_t i = 0; i < this->N; ++i) {
                this->v[i][w + 1] = 0;
            }

            this->increase_counter();
        } else {
            for (std::size_t i = 0; i < this->N; ++i) {
                this->v[i][w] = this->v[i][w - 1];
            }
            std::fill(this->Y.begin(), this->Y.end(), 0);
        }
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            auto const result = this->Y[this->j++];
            if (this->j == this->N) {
                this->generate();
                this->increase_counter();
                this->j = 0;
            }

            return result;
        } else {
            if (this->j == 0) {
                this->j = this->N;
                this->decrease_counter();
                this->generate();
            }
            --this->j;
            return this->Y[this->j];
        }
    };

    inline void discard(unsigned long long z) {

        if constexpr (skip_first) {
            this->decrease_counter();
        }

        constexpr auto digits = std::numeric_limits<UIntType>::digits;
        using upgraded_type =
            std::conditional_t<digits <= 32, std::uint64_t, uint128>;

        upgraded_type zup{z};

        upgraded_type const ref_n_sims =
            (static_cast<upgraded_type>(this->m_seq_num) *
             static_cast<upgraded_type>(this->N)) +
            (static_cast<upgraded_type>(this->j));

        upgraded_type const full_cycle =
            (static_cast<upgraded_type>(sobol_engine::max()) +
             upgraded_type(1)) *
            static_cast<upgraded_type>(this->N);
        upgraded_type z_rem = (zup % full_cycle);

        z_rem += ref_n_sims;
        // we might have gone over when adding curent paths
        z_rem %= full_cycle;

        upgraded_type Nup{this->N};

        using adhoc::lldiv;
        using std::lldiv;
        auto const res = lldiv(z_rem, Nup);

        this->m_seq_num = static_cast<UIntType>(res.quot);

        // double cast because for some reason uint128 cannot have conversion
        // to both std::size_t and std::uint64_t
        this->j = static_cast<std::size_t>(static_cast<std::uint64_t>(res.rem));

        if constexpr (skip_first) {
            this->increase_counter();
        }

        // generate from sequence number
        if constexpr (skip_first) {
            if (this->m_seq_num == 0) {
                for (std::size_t i = 0; i < this->N; ++i) {
                    this->Y[i] = this->v[i][0];
                }
                return;
            }
        }

        std::fill(this->Y.begin(), this->Y.end(), 0);

        // skip to sequence number seq_num
        for (std::size_t i = 0; i < this->N; ++i) {
            // figure out Gray code of seq_num
            UIntType seq_num_gray = this->m_seq_num ^ (this->m_seq_num >> 1U);
            // construct x[i] based on seq_num_gray

            for (std::size_t ii = 0; ii < w; ++ii) {
                if (seq_num_gray == 0) {
                    break;
                }

                if ((seq_num_gray & 1U) != 0U) {
                    this->Y[i] ^= this->v[i][ii];
                }

                seq_num_gray >>= 1U;
            }
        }
    }

    static constexpr auto min() -> UIntType {
        if constexpr (skip_first) {
            return static_cast<UIntType>(1U);
        } else {
            return static_cast<UIntType>(0U);
        }
    }

    static constexpr auto max() -> UIntType { return mask<UIntType, w>(); }

    auto operator==(const sobol_engine &rhs) const -> bool {
        return (this->j == rhs.j) && (this->N == rhs.N) &&
               (this->m_seq_num == rhs.m_seq_num) && (this->v == rhs.v) &&
               (this->Y == rhs.Y);
    }

  private:
    std::size_t j{0};
    std::size_t N{1};

    UIntType m_seq_num{static_cast<UIntType>(0)}; // sequence number

    std::vector<
        std::array<UIntType, w + 1 + static_cast<std::size_t>(skip_first)>>
        v; // direction numbers
    std::vector<UIntType> Y;
};

template <bool skip_first = true>
using sobol = sobol_engine<std::uint_fast32_t, 32, skip_first>;

template <bool skip_first = true>
using sobol64 = sobol_engine<std::uint_fast64_t, 64, skip_first>;

} // namespace adhoc

#endif // BRNG_SOBOL_ENGINE
