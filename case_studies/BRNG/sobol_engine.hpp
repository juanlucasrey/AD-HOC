#ifndef CASE_STUDIES_BRNG_SOBOL_ENGINE
#define CASE_STUDIES_BRNG_SOBOL_ENGINE

#include "bit_find.hpp"
#include "sobol_tables/joe_kuo_old_1111.hpp"
#include "sobol_tables/joe_kuo_other_0_7600.hpp"
#include "sobol_tables/joe_kuo_other_2_3900.hpp"
#include "sobol_tables/joe_kuo_other_3_7300.hpp"
#include "sobol_tables/joe_kuo_other_4_5600.hpp"
#include "sobol_tables/new_joe_kuo_5_21201.hpp"
#include "sobol_tables/new_joe_kuo_6_21201.hpp"
#include "sobol_tables/new_joe_kuo_7_21201.hpp"
#include "uint128.hpp"

#include <array>
#include <bitset>
#include <cstdint>
#include <exception>
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

constexpr auto poly(qrng_table table, std::size_t i) -> unsigned int {
    switch (table) {
    case qrng_table::joe_kuo_old_1111:
        return joe_kuo_old_1111_poly[i];
    case qrng_table::joe_kuo_other_0_7600:
        return joe_kuo_other_0_7600_poly[i];
    case qrng_table::joe_kuo_other_2_3900:
        return joe_kuo_other_2_3900_poly[i];
    case qrng_table::joe_kuo_other_3_7300:
        return joe_kuo_other_3_7300_poly[i];
    case qrng_table::joe_kuo_other_4_5600:
        return joe_kuo_other_4_5600_poly[i];
    case qrng_table::new_joe_kuo_5_21201:
        return new_joe_kuo_5_21201_poly[i];
    case qrng_table::new_joe_kuo_6_21201:
        return new_joe_kuo_6_21201_poly[i];
    case qrng_table::new_joe_kuo_7_21201:
        return new_joe_kuo_7_21201_poly[i];
    }
    return 0;
}

template <qrng_table t>
auto v_init(std::size_t i, std::size_t j) -> unsigned int;

template <>
inline auto v_init<qrng_table::joe_kuo_old_1111>(std::size_t i, std::size_t j)
    -> unsigned int {
    return static_cast<unsigned int>(joe_kuo_old_1111_init[i][j]);
}

template <>
inline auto v_init<qrng_table::joe_kuo_other_0_7600>(std::size_t i,
                                                     std::size_t j)
    -> unsigned int {
    if (j < joe_kuo_other_0_7600_cut_degree) {
        return static_cast<unsigned int>(joe_kuo_other_0_7600_init[i][j]);
    }

    if (i < joe_kuo_other_0_7600_cut_dimension) {
        return 0;
    }

    return joe_kuo_other_0_7600_init_last[i -
                                          joe_kuo_other_0_7600_cut_dimension]
                                         [j - joe_kuo_other_0_7600_cut_degree];
}

template <>
inline auto v_init<qrng_table::joe_kuo_other_2_3900>(std::size_t i,
                                                     std::size_t j)
    -> unsigned int {
    return static_cast<unsigned int>(joe_kuo_other_2_3900_init[i][j]);
}

template <>
inline auto v_init<qrng_table::joe_kuo_other_3_7300>(std::size_t i,
                                                     std::size_t j)
    -> unsigned int {
    if (j < joe_kuo_other_3_7300_cut_degree) {
        return static_cast<unsigned int>(joe_kuo_other_3_7300_init[i][j]);
    }

    if (i < joe_kuo_other_3_7300_cut_dimension) {
        return 0;
    }

    return joe_kuo_other_3_7300_init_last[i -
                                          joe_kuo_other_3_7300_cut_dimension]
                                         [j - joe_kuo_other_3_7300_cut_degree];
}

template <>
inline auto v_init<qrng_table::joe_kuo_other_4_5600>(std::size_t i,
                                                     std::size_t j)
    -> unsigned int {
    return static_cast<unsigned int>(joe_kuo_other_4_5600_init[i][j]);
}

template <>
inline auto v_init<qrng_table::new_joe_kuo_5_21201>(std::size_t i,
                                                    std::size_t j)
    -> unsigned int {
    if (j < new_joe_kuo_5_21201_cut_degree) {
        return static_cast<unsigned int>(new_joe_kuo_5_21201_init[i][j]);
    }

    if (i < new_joe_kuo_5_21201_cut_dimension) {
        return 0;
    }

    return new_joe_kuo_5_21201_init_last[i - new_joe_kuo_5_21201_cut_dimension]
                                        [j - new_joe_kuo_5_21201_cut_degree];
}

template <>
inline auto v_init<qrng_table::new_joe_kuo_6_21201>(std::size_t i,
                                                    std::size_t j)
    -> unsigned int {
    if (j < new_joe_kuo_6_21201_cut_degree) {
        return static_cast<unsigned int>(new_joe_kuo_6_21201_init[i][j]);
    }

    if (i < new_joe_kuo_6_21201_cut_dimension) {
        return 0;
    }

    return new_joe_kuo_6_21201_init_last[i - new_joe_kuo_6_21201_cut_dimension]
                                        [j - new_joe_kuo_6_21201_cut_degree];
}

template <>
inline auto v_init<qrng_table::new_joe_kuo_7_21201>(std::size_t i,
                                                    std::size_t j)
    -> unsigned int {
    if (j < new_joe_kuo_7_21201_cut_degree) {
        return static_cast<unsigned int>(new_joe_kuo_7_21201_init[i][j]);
    }

    if (i < new_joe_kuo_7_21201_cut_dimension) {
        return 0;
    }

    return new_joe_kuo_7_21201_init_last[i - new_joe_kuo_7_21201_cut_dimension]
                                        [j - new_joe_kuo_7_21201_cut_degree];
}

} // namespace detail

template <class UIntType, std::size_t w, bool skip_first = true,
          qrng_table t = qrng_table::joe_kuo_old_1111>
class sobol_engine final {
    static_assert(w <= std::numeric_limits<UIntType>::digits);
    static_assert(w <= 64);

    inline void increase_counter() {
        constexpr auto in_mask = sobol_engine::max();
        this->m_seq_num = (this->m_seq_num + 1) & in_mask;
    }

    inline void decrease_counter() {
        constexpr auto in_mask = sobol_engine::max();
        this->m_seq_num =
            this->m_seq_num == 0 ? in_mask : (this->m_seq_num - 1);
    }

    inline void generate() {

        constexpr auto digits = std::numeric_limits<UIntType>::digits;
        using ceiling_type =
            std::conditional_t<digits <= 32, std::uint32_t, std::uint64_t>;

        auto const c = findLeastSignificantZero(
            static_cast<ceiling_type>(this->m_seq_num));
        for (std::size_t i = 0; i < this->N; ++i) {
            this->Y[i] ^= this->v[i][c];
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

        for (std::size_t i = 1; i < this->N; ++i) {
            unsigned int poly_e = detail::poly(t, i);
            std::size_t const poly_d = highestBitIndex(poly_e);
            std::bitset<w> includ(poly_e);

            std::size_t ii = 0;
            for (; ii < poly_d; ++ii) {
                this->v[i][ii] = detail::v_init<t>(i, ii);
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

        for (std::size_t i = 0; i < this->N; ++i) {
            this->v[i][w] = this->v[i][w - 1];
        }

        this->generate_from_num();
    }

    inline auto operator()() -> result_type {
        auto const result = this->Y[this->j++];
        if (this->j == this->N) {
            this->generate();
            this->increase_counter();
            this->j = 0;
        }

        return result;
    };

    inline void generate_from_num() {
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
        this->j = static_cast<std::size_t>(res.rem);

        if constexpr (skip_first) {
            this->increase_counter();
        }

        this->generate_from_num();
    }

    static constexpr auto min() -> UIntType {
        return static_cast<UIntType>(0U);
    }

    static constexpr auto max() -> UIntType {
        constexpr UIntType ach = static_cast<UIntType>(1U)
                                 << static_cast<std::size_t>(w - 1U);

        return ach | (ach - static_cast<UIntType>(1U));

        // equivalent to the following but shift does not overflow
        // return (static_cast<result_type>(1U) << w) - 1U;
    }

    auto operator==(const sobol_engine &rhs) const -> bool {
        return (this->j == rhs.j) && (this->N == rhs.N) &&
               (this->m_seq_num == rhs.m_seq_num) && (this->v == rhs.v) &&
               (this->Y == rhs.Y);
    }

  private:
    std::size_t j{0};
    std::size_t N{1};

    UIntType m_seq_num{
        static_cast<unsigned int>(skip_first)}; // sequence number

    std::vector<std::array<UIntType, w + 1>> v; // direction numbers
    std::vector<UIntType> Y;
};

using sobol = sobol_engine<std::uint_fast32_t, 32>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_SOBOL_ENGINE
