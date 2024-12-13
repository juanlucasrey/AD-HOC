#ifndef CASE_STUDIES_BRNG_PHILOX_ENGINE
#define CASE_STUDIES_BRNG_PHILOX_ENGINE

#include "uint128.hpp"

#include <array>
#include <cstdint>
#include <limits>

namespace adhoc {

template <typename UIntType, std::size_t w, std::size_t n, std::size_t r,
          UIntType... consts>
class philox_engine {
  private:
    static_assert(sizeof...(consts) == n);
    static_assert(n == 2 || n == 4 || n == 8 || n == 16);
    static_assert(0 < r);
    static_assert(0 < w && w <= std::numeric_limits<UIntType>::digits);

    inline void increase_counter() {
        constexpr auto in_mask = max();
        std::size_t i = 0;
        do {
            this->X[i] = (this->X[i] + 1) & in_mask;
            ++i;
        } while (i < n && !this->X[i - 1]);
    }

    inline void decrease_counter() {
        constexpr auto in_mask = max();
        std::size_t i = 0;
        do {
            this->X[i] = (this->X[i] - 1) & in_mask;
            ++i;
        } while (i < n && (this->X[i - 1] == in_mask));
    }

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t word_count = n;
    static constexpr std::size_t round_count = r;
    static constexpr result_type default_seed = 20111115U;

    static constexpr auto min() -> result_type {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> result_type {
        return (result_type(~(result_type(0))) >>
                (std::numeric_limits<result_type>::digits - w));

        // equivalent to the following but shift does not overflow
        // return (static_cast<result_type>(1U) << w) - 1U;
    }

    philox_engine() : philox_engine(default_seed) {}

    explicit philox_engine(result_type value) {
        constexpr result_type mask = max();
        this->K[0] = value & mask;
        this->decrease_counter();
        this->generate<2>();
        this->increase_counter();
    }

    template <class SeedSeq> explicit philox_engine(SeedSeq &seq) {
        constexpr std::size_t p = (w - 1) / 32 + 1;
        constexpr std::size_t n_half = n / 2;

        std::array<std::uint_least32_t, n_half * p> a;
        seq.generate(a.begin(), a.end());

        constexpr result_type mask = max();
        auto iter = a.begin();
        for (std::size_t n_idx = 0; n_idx < n_half; n_idx++) {
            UIntType val = 0;
            for (std::size_t p_idx = 0; p_idx < p; ++p_idx) {
                val += static_cast<UIntType>(*iter++) << 32 * p_idx;
            }
            this->K[n_idx] = val & mask;
        }
    }

    void set_counter(const std::array<result_type, n> &c) {
        constexpr result_type mask = max();
        for (std::size_t i = 0; i < n; ++i) {
            this->X[i] = c[i] & mask;
        }
        this->j = n - 1;
        this->decrease_counter();
        this->generate<2>();
        this->increase_counter();
    }
    // version 0: STD
    // version 1: correction, working
    // version 2: correction, efficient
    template <bool FwdDirection = true, std::size_t Version = 2>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            ++this->j;
            if (this->j == n) {
                this->generate<Version>();
                this->increase_counter();
                this->j = 0;
            }
            return Y[this->j];
        } else {
            const result_type result = Y[this->j];
            if (this->j == 0) {
                this->j = n;
                this->decrease_counter();
                this->decrease_counter();
                this->generate<Version>();
                this->increase_counter();
            }
            --this->j;
            return result;
        }
    }

    // version 0: STD
    // version 1: correction, working
    // version 2: correction, efficient
    template <std::size_t Version = 2>
    inline void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()<true, Version>();
        }
    }

    auto operator==(const philox_engine &rhs) const -> bool {
        return (this->X == rhs.X) && (this->K == rhs.K) && (this->Y == rhs.Y) &&
               (this->j == rhs.j);
    }

  private:
    template <std::unsigned_integral U>
    inline auto mulhilo(U a, U b) -> std::pair<U, U> {

        using upgraded_type = std::conditional_t<
            w <= 8, std::uint_fast16_t,
            std::conditional_t<
                w <= 16, std::uint_fast32_t,
                std::conditional_t<w <= 32, std::uint_fast64_t, uint128>>>;

        const upgraded_type ab =
            static_cast<upgraded_type>(a) * static_cast<upgraded_type>(b);
        return {static_cast<U>(ab >> w), static_cast<U>(ab) & this->max()};
    }

    template <std::size_t Version> inline void generate() {
        constexpr auto in_mask = max();
        constexpr std::array<UIntType, n> consts_arr{consts...};
        if constexpr (n == 2) {

            if constexpr (Version == 0) {
                // version 0: STD
                result_type S0 = this->X[0];
                result_type S1 = this->X[1];
                result_type K0 = this->K[0];
                for (size_t i = 0; i < round_count; ++i) {
                    auto [hi, lo] = this->mulhilo(S1, consts_arr[0]);
                    S0 = lo;
                    S1 = hi ^ K0 ^ S0;
                    K0 = (K0 + consts_arr[1]) & in_mask;
                }
                this->Y[0] = S0;
                this->Y[1] = S1;
            } else {
                // version 1: correction, working
                // version 2: correction, efficient
                result_type S0 = this->X[0];
                result_type S1 = this->X[1];
                result_type K0 = this->K[0];
                for (size_t i = 0; i < round_count; ++i) {
                    auto [hi, lo] = this->mulhilo(S0, consts_arr[0]);
                    S0 = hi ^ K0 ^ S1;
                    S1 = lo;
                    K0 = (K0 + consts_arr[1]) & in_mask;
                }
                this->Y[0] = S0;
                this->Y[1] = S1;
            }

        } else if constexpr (n == 4) {
            if constexpr (Version == 0) {
                // version 0: STD
                result_type S0 = this->X[0];
                result_type S1 = this->X[1];
                result_type S2 = this->X[2];
                result_type S3 = this->X[3];
                result_type K0 = this->K[0];
                result_type K1 = this->K[1];
                for (size_t i = 0; i < r; ++i) {
                    result_type V0 = S0;
                    result_type V1 = S3;
                    result_type V2 = S2;
                    result_type V3 = S1;

                    auto [hi1, lo1] = this->mulhilo(V1, consts_arr[0]);
                    auto [hi3, lo3] = this->mulhilo(V3, consts_arr[2]);
                    S0 = lo1;
                    S1 = hi1 ^ K0 ^ V0;
                    S2 = lo3;
                    S3 = hi3 ^ K1 ^ V2;
                    K0 = (K0 + consts_arr[1]) & in_mask;
                    K1 = (K1 + consts_arr[3]) & in_mask;
                }
                this->Y[0] = S0;
                this->Y[1] = S1;
                this->Y[2] = S2;
                this->Y[3] = S3;
            } else if constexpr (Version == 1) {
                // version 1: correction, working
                result_type S0 = this->X[0];
                result_type S1 = this->X[1];
                result_type S2 = this->X[2];
                result_type S3 = this->X[3];
                result_type K0 = this->K[0];
                result_type K1 = this->K[1];
                for (size_t i = 0; i < r; ++i) {
                    result_type V0 = S2;
                    result_type V1 = S1;
                    result_type V2 = S0;
                    result_type V3 = S3;

                    auto [hi0, lo0] = this->mulhilo(
                        V0, consts_arr[2]); // multiplier inverted!
                    auto [hi2, lo2] = this->mulhilo(
                        V2, consts_arr[0]); // multiplier inverted!

                    S0 = hi0 ^ K0 ^ V1;
                    S1 = lo0;
                    S2 = hi2 ^ K1 ^ V3;
                    S3 = lo2;
                    K0 = (K0 + consts_arr[1]) & in_mask;
                    K1 = (K1 + consts_arr[3]) & in_mask;
                }
                this->Y[0] = S0;
                this->Y[1] = S1;
                this->Y[2] = S2;
                this->Y[3] = S3;
            } else {
                // version 2: correction, efficient
                result_type S0 = this->X[0];
                result_type S1 = this->X[1];
                result_type S2 = this->X[2];
                result_type S3 = this->X[3];
                result_type K0 = this->K[0];
                result_type K1 = this->K[1];
                for (size_t i = 0; i < r; ++i) {
                    auto [hi0, lo0] = this->mulhilo(
                        S2, consts_arr[2]); // multiplier inverted!
                    auto [hi2, lo2] = this->mulhilo(
                        S0, consts_arr[0]); // multiplier inverted!
                    S0 = hi0 ^ K0 ^ S1;
                    S1 = lo0;
                    S2 = hi2 ^ K1 ^ S3;
                    S3 = lo2;
                    K0 = (K0 + consts_arr[1]) & in_mask;
                    K1 = (K1 + consts_arr[3]) & in_mask;
                }
                this->Y[0] = S0;
                this->Y[1] = S1;
                this->Y[2] = S2;
                this->Y[3] = S3;
            }
        }
    }

    std::array<UIntType, n> X{0};
    std::array<UIntType, n / 2> K{0};
    std::array<UIntType, n> Y{0};
    std::size_t j{n - 1};
};

using philox4x32 = philox_engine<std::uint_fast32_t, 32, 4, 10, 0xD2511F53,
                                 0x9E3779B9, 0xCD9E8D57, 0xBB67AE85>;

using philox4x64 =
    philox_engine<std::uint_fast64_t, 64, 4, 10, 0xD2E7470EE14C6C93,
                  0x9E3779B97F4A7C15, 0xCA5A826395121157, 0xBB67AE8584CAA73B>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_PHILOX_ENGINE
