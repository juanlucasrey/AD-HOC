#ifndef CASE_STUDIES_BRNG_MERSENNE_TWISTER_ENGINE
#define CASE_STUDIES_BRNG_MERSENNE_TWISTER_ENGINE

#include <array>
#include <cstdint>
#include <limits>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t n, std::size_t m,
          std::size_t r, UIntType a, std::size_t u, UIntType d, std::size_t s,
          UIntType b, std::size_t t, UIntType c, std::size_t l, UIntType f>
class mersenne_twister_engine final {
  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t state_size = n;
    static constexpr std::size_t shift_size = m;
    static constexpr std::size_t mask_bits = r;
    static constexpr UIntType xor_mask = a;
    static constexpr std::size_t tempering_u = u;
    static constexpr UIntType tempering_d = d;
    static constexpr std::size_t tempering_s = s;
    static constexpr UIntType tempering_b = b;
    static constexpr std::size_t tempering_t = t;
    static constexpr UIntType tempering_c = c;
    static constexpr std::size_t tempering_l = l;
    static constexpr UIntType initialization_multiplier = f;
    static constexpr UIntType default_seed = 5489U;

    explicit mersenne_twister_engine(result_type value = default_seed) {
        this->i = 0;
        this->x[this->i++] = value;

        for (this->i = 1; this->i < n; this->i++) {
            this->x[this->i] = (f * (this->x[this->i - 1] ^
                                     (this->x[this->i - 1] >> (w - 2))) +
                                static_cast<UIntType>(this->i));
        }

        // this short step ensures consistency when the full_cycle is completed,
        // making simul.forward[full_cycle - 2] = simul.backward[1]
        constexpr UIntType upper_mask = (~static_cast<UIntType>(0)) << r;
        constexpr UIntType lower_mask = ~upper_mask;

        UIntType tmp;
        UIntType y;
        y = this->x[n - 1] ^ this->x[m - 1];
        tmp = y >> (std::numeric_limits<UIntType>::digits - 1);
        y = ((y ^ (tmp * a)) << 1) | tmp;
        this->x[0] = (this->x[0] & upper_mask) | (y & lower_mask);
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        UIntType z;

        constexpr UIntType upper_mask = (~static_cast<UIntType>(0)) << r;
        constexpr UIntType lower_mask = ~upper_mask;
        constexpr std::size_t n_minus_m = n - m;

        if constexpr (FwdDirection) {
            if (this->i == n) {
                UIntType y;
                std::size_t kk;

                for (kk = 0; kk < n_minus_m; kk++) {
                    y = (this->x[kk] & upper_mask) |
                        (this->x[kk + 1] & lower_mask);
                    this->x[kk] = this->x[kk + m] ^ (y >> 1) ^
                                  ((this->x[kk + 1] & 1U) * a);
                }

                for (; kk < n - 1; kk++) {
                    y = (this->x[kk] & upper_mask) |
                        (this->x[kk + 1] & lower_mask);
                    this->x[kk] = this->x[kk - n_minus_m] ^ (y >> 1) ^
                                  ((this->x[kk + 1] & 1U) * a);
                }

                y = (this->x[n - 1] & upper_mask) | (this->x[0] & lower_mask);
                this->x[n - 1] =
                    this->x[m - 1] ^ (y >> 1) ^ ((this->x[0] & 1U) * a);
                this->i = 0;
            }
            z = this->x[this->i++];
        } else {
            z = this->x[(--this->i)];
            if (this->i == 0) {
                UIntType tmp;
                UIntType y;
                std::size_t kk;

                y = this->x[n - 1] ^ this->x[m - 1];
                tmp = y >> (std::numeric_limits<UIntType>::digits - 1);
                y = ((y ^ (tmp * a)) << 1) | tmp;
                this->x[n - 1] =
                    (y & upper_mask) | (this->x[n - 1] & lower_mask);
                this->x[0] = (this->x[0] & upper_mask) | (y & lower_mask);

                for (kk = n - 1; kk-- > n_minus_m;) {
                    y = this->x[kk] ^ this->x[kk - n_minus_m];
                    tmp = y >> (std::numeric_limits<UIntType>::digits - 1);
                    y = ((y ^ (tmp * a)) << 1) | tmp;
                    this->x[kk] = (y & upper_mask) | (this->x[kk] & lower_mask);
                    this->x[kk + 1] =
                        (this->x[kk + 1] & upper_mask) | (y & lower_mask);
                }

                for (kk = n_minus_m; kk-- > 0;) {
                    y = this->x[kk] ^ this->x[kk + m];
                    tmp = y >> (std::numeric_limits<UIntType>::digits - 1);
                    y = ((y ^ (tmp * a)) << 1) | tmp;
                    this->x[kk] = (y & upper_mask) | (this->x[kk] & lower_mask);
                    this->x[kk + 1] =
                        (this->x[kk + 1] & upper_mask) | (y & lower_mask);
                }

                y = this->x[n - 1] ^ this->x[m - 1];
                tmp = y >> (std::numeric_limits<UIntType>::digits - 1);
                y = ((y ^ (tmp * a)) << 1) | tmp;
                this->x[0] = (this->x[0] & upper_mask) | (y & lower_mask);
                this->i = n;
            }
        }

        z ^= ((z >> u) & d);
        z ^= ((z << s) & b);
        z ^= ((z << t) & c);
        z ^= (z >> l);
        return z;
    }

    static constexpr auto min() -> UIntType {
        return std::numeric_limits<UIntType>::min();
    }
    static constexpr auto max() -> UIntType {
        return std::numeric_limits<UIntType>::max();
    }

  private:
    std::array<UIntType, n> x{0};
    std::size_t i{0};
};

// using mt19937 =
//     mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31, 0x9908b0df,
//                             11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000,
//                             18, 1812433253>;
using mt19937 =
    mersenne_twister_engine<std::uint32_t, 32, 624, 397, 31, 0x9908b0df, 11,
                            0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18,
                            1812433253>;

using mt19937_64 =
    mersenne_twister_engine<std::uint_fast64_t, 64, 312, 156, 31,
                            0xb5026f5aa96619e9, 29, 0x5555555555555555, 17,
                            0x71d67fffeda60000, 37, 0xfff7eee000000000, 43,
                            6364136223846793005>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_MERSENNE_TWISTER_ENGINE
