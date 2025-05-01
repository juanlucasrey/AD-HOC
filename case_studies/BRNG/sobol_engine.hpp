#ifndef CASE_STUDIES_BRNG_SOBOL_ENGINE
#define CASE_STUDIES_BRNG_SOBOL_ENGINE

#include "sobol_tables/joe_kuo_old_1111.hpp"
#include "sobol_tables/joe_kuo_other_0_7600.hpp"
#include "sobol_tables/joe_kuo_other_2_3900.hpp"
#include "sobol_tables/joe_kuo_other_3_7300.hpp"
#include "sobol_tables/joe_kuo_other_4_5600.hpp"
#include "sobol_tables/new_joe_kuo_5_21201.hpp"
#include "sobol_tables/new_joe_kuo_6_21201.hpp"
#include "sobol_tables/new_joe_kuo_7_21201.hpp"

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

inline auto keepHighestBit(unsigned int n) -> unsigned int {
    n |= (n >> 1);
    n |= (n >> 2);
    n |= (n >> 4);
    n |= (n >> 8);
    n |= (n >> 16);
    return n - (n >> 1);
}

inline auto highestBitIndex(unsigned int b) -> unsigned int {
    constexpr unsigned int deBruijnMagic = 0x06EB14F9;
    constexpr std::array<unsigned int, 32> deBruijnTable = {
        0,  1,  16, 2,  29, 17, 3,  22, 30, 20, 18, 11, 13, 4, 7,  23,
        31, 15, 28, 21, 19, 10, 12, 6,  14, 27, 9,  5,  26, 8, 25, 24,
    };
    return deBruijnTable[(keepHighestBit(b) * deBruijnMagic) >> 27];
}

} // namespace detail

template <class UIntType, std::size_t w, bool skip_first = true,
          qrng_table t = qrng_table::joe_kuo_old_1111>
class sobol_engine final {
    static_assert(w <= std::numeric_limits<UIntType>::digits);

    inline void increase_counter() {
        constexpr auto in_mask = sobol_engine::max();
        this->m_seq_num = (this->m_seq_num + 1) & in_mask;
    }

    inline void decrease_counter() {
        constexpr auto in_mask = sobol_engine::max();
        this->m_seq_num =
            this->m_seq_num == 0 ? in_mask : (this->m_seq_num - 1);
    }

    inline auto findLeastSignificantZeroNaive(UIntType n) {
        UIntType c = 0;
        while ((n & 1)) {
            n >>= 1;
            ++c;
        }
        return c;
    }

    // Returns the position of the least significant 0 bit (0-based index from
    // right) Returns 32 if no 0 bits exist (i.e., input is all 1s)
    inline auto findLeastSignificantZero32(uint32_t n) -> uint32_t {
        // Invert bits: 1s become 0s and 0s become 1s
        n = ~n;

        static const int deBruijnTable[32] = {
            0,  1,  28, 2,  29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4,  8,
            31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6,  11, 5,  10, 9};

        uint32_t isolated = n & -n;
        // Prevent multiplication by 0 and ensure valid table index
        uint32_t adjusted = isolated | (isolated == 0);
        uint32_t index = ((adjusted * 0x077CB531U) >> 27);
        // Use bitwise ops to select between table value and 32
        uint32_t is_zero = (isolated == 0);
        return (deBruijnTable[index] & ~(-is_zero)) | (32 & (-is_zero));
    }

    // Returns the position of the least significant 0 bit (0-based index from
    // right) Returns 64 if no 0 bits exist (i.e., input is all 1s)
    inline auto findLeastSignificantZero64(uint64_t n) -> uint64_t {
        // Invert bits: 1s become 0s and 0s become 1s
        n = ~n;

        static const int deBruijnTable[64] = {
            0,  1,  2,  7,  3,  13, 8,  19, 4,  25, 14, 33, 9,  36, 20, 44,
            5,  28, 26, 38, 15, 46, 34, 51, 10, 39, 37, 52, 21, 54, 45, 57,
            63, 6,  12, 18, 24, 32, 35, 43, 27, 47, 50, 56, 62, 11, 17, 23,
            31, 42, 49, 55, 61, 16, 22, 30, 41, 48, 60, 29, 40, 59, 58, 53};

        uint64_t isolated = n & -n;
        // Prevent multiplication by 0 and ensure valid table index
        uint64_t adjusted = isolated | (isolated == 0);
        uint64_t index = ((adjusted * 0x03F79D71B4CB0A89ULL) >> 58);
        // Use bitwise ops to select between table value and 64
        uint64_t is_zero = (isolated == 0);
        return (deBruijnTable[index] & ~(-is_zero)) | (64 & (-is_zero));
    }

    inline void generate() {
        auto const c =
            findLeastSignificantZero32(static_cast<uint32_t>(this->m_seq_num));
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

        for (unsigned int i = 1; i < this->N; ++i) {
            unsigned int poly_e = detail::poly(t, i);
            unsigned int const poly_d = detail::highestBitIndex(poly_e);
            std::bitset<w> includ(poly_e);

            unsigned int ii = 0;
            for (; ii < poly_d; ++ii) {
                this->v[i][ii] = detail::v_init<t>(i, ii);
            }

            for (; ii < w; ++ii) {
                auto &v_new = this->v[i][ii];
                v_new = this->v[i][ii - poly_d];
                // unsigned int v_new = this->v[i][ii - poly_d];
                unsigned int kk = 1;

                for (unsigned int k = 0; k < poly_d; ++k) {
                    kk *= 2;
                    if (includ[poly_d - k - 1]) {
                        v_new ^= (kk * v[i][ii - k - 1]);
                    }
                }
            }
        }

        unsigned int kk = 1;
        for (unsigned int i = w - 1; i-- > 0;) {
            kk *= 2;

            for (unsigned int ii = 0; ii < this->N; ++ii) {
                this->v[ii][i] *= kk;
            }
        }

        for (unsigned int i = 0; i < this->N; ++i) {
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
        for (unsigned int i = 0; i < this->N; ++i) {
            // figure out Gray code of seq_num
            unsigned int seq_num_gray =
                this->m_seq_num ^ (this->m_seq_num >> 1U);
            // construct x[i] based on seq_num_gray

            for (unsigned int ii = 0; ii < w; ++ii) {
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

        uint64_t const ref_n_sims = (static_cast<uint64_t>(this->m_seq_num) *
                                     static_cast<uint64_t>(this->N)) +
                                    (static_cast<uint64_t>(this->j));
        uint64_t const full_cycle =
            static_cast<uint64_t>(sobol_engine::max() + 1) *
            static_cast<uint64_t>(this->N);
        uint64_t z_rem = (z % full_cycle);

        z_rem += ref_n_sims;
        // we might have gone over when adding curent paths
        z_rem %= full_cycle;

        this->m_seq_num = static_cast<unsigned int>(z_rem / this->N);
        this->j =
            static_cast<unsigned int>(z_rem % static_cast<uint64_t>(this->N));

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
    unsigned int j{0};
    unsigned int N{1};

    unsigned int m_seq_num{
        static_cast<unsigned int>(skip_first)}; // sequence number

    std::vector<std::array<unsigned int, w + 1>> v; // direction numbers
    std::vector<UIntType> Y;
};

using sobol = sobol_engine<std::uint_fast32_t, 32>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_SOBOL_ENGINE
