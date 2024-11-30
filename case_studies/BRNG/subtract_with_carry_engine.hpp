#ifndef CASE_STUDIES_BRNG_SUBTRACT_WITH_CARRY_ENGINE
#define CASE_STUDIES_BRNG_SUBTRACT_WITH_CARRY_ENGINE

#include "linear_congruential_engine.hpp"

#include "../../include/combinatorics/pow.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include <iostream>

namespace adhoc {

// G. Marsaglia, A. Zaman, A new class of random number generators, The Annals
// of Applied Probability 1 (3) (1991) 462–480. URL
// http://www.jstor.org/stable/2959748

// on page 11 of the article (page 472 of the publication) it says: "Finding
// conditions which ensure that a seed vector is periodic is pretty much an
// acedemic exercise, done out of curiosity"

// the corrected=true template parameter ensures that our RNG is strictly
// periodic.

// Marsaglia et al. pointed at this when he wrote on page 3 (page 464 of the
// publication) "since f has an inverse [...] the sequence is strictly periodic"
// corrected = true ensures that f has an inverse and is therefore strictly
// periodic.

template <class UIntType, std::size_t w, std::size_t s, std::size_t r,
          bool corrected = false>
class subtract_with_carry_engine final {
  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t short_lag = s;
    static constexpr std::size_t long_lag = r;
    static constexpr std::uint_least32_t default_seed = 19780503U;

    static constexpr result_type modulus = static_cast<result_type>(1U) << w;

    explicit subtract_with_carry_engine(result_type value = default_seed) {
        linear_congruential_engine<std::uint_least32_t, 40014U, 0U, 2147483563U>
            e(value == 0U ? default_seed : value);
        constexpr auto n = std::size_t(w / 32) + 1;

        for (std::size_t j = 0; j < r; j++) {
            UIntType val = 0;
            for (std::size_t k = 0; k < n; ++k) {
                val += static_cast<UIntType>(e()) << 32 * k;
            }
            x[j] = val % modulus;
        }

        k = 0;

        if constexpr (corrected) {
            std::size_t k_prev = k == 0 ? (long_lag - 1) : (k - 1);
            std::size_t short_index_prev = (k_prev < short_lag)
                                               ? (k_prev + long_lag - short_lag)
                                               : (k_prev - short_lag);

            UIntType temp_prev = x[k_prev];
            if (temp_prev > x[short_index_prev]) {
                temp_prev = modulus - temp_prev + x[short_index_prev];
            } else {
                temp_prev = x[short_index_prev] - temp_prev;
            }

            if (x[short_index_prev] >= temp_prev) {
                carry = 0;
            } else {
                carry = 1;
            }
        } else {
            carry = (x[long_lag - 1] == 0);
        }
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        if constexpr (FwdDirection) {
            std::size_t short_index =
                (k < short_lag) ? (k + long_lag - short_lag) : (k - short_lag);

            UIntType temp = x[k] + carry;
            if (x[short_index] >= temp) {
                x[k] = x[short_index] - temp;
                carry = 0;
            } else {
                x[k] = modulus - temp + x[short_index];
                carry = 1;
            }

            UIntType result = x[k];
            ++k;
            if (k >= long_lag) {
                k = 0;
            }
            return result;
        } else {
            if (k == 0) {
                k = long_lag - 1;
            } else {
                --k;
            }
            UIntType result = x[k];

            std::size_t short_index =
                (k < short_lag) ? (k + long_lag - short_lag) : (k - short_lag);

            UIntType temp;
            if (carry) {
                temp = modulus - x[k] + x[short_index];
            } else {
                temp = x[short_index] - x[k];
            }

            std::size_t k_prev = k == 0 ? (long_lag - 1) : (k - 1);
            std::size_t short_index_prev = (k_prev < short_lag)
                                               ? (k_prev + long_lag - short_lag)
                                               : (k_prev - short_lag);

            UIntType temp_prev = x[k_prev];
            if (temp_prev > x[short_index_prev]) {
                temp_prev = modulus - temp_prev + x[short_index_prev];
            } else {
                temp_prev = x[short_index_prev] - temp_prev;
            }

            if (x[short_index_prev] >= temp_prev) {
                carry = 0;
            } else {
                carry = 1;
            }

            x[k] = temp - carry;

            return result;
        }
    }

    void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }

    static constexpr auto period() -> unsigned long long {
        // return std::pow(w, r) - std::pow(w, s);
        return adhoc::detail::pow<r>(w) - adhoc::detail::pow<s>(w);
    }

    static constexpr auto min() -> UIntType {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> UIntType {
        return (static_cast<result_type>(1U) << w) - 1;
    }

    auto operator==(const subtract_with_carry_engine &rhs) const -> bool {
        return (this->carry == rhs.carry) && (this->k == rhs.k) &&
               (this->x == rhs.x);
    }

  private:
    std::array<UIntType, long_lag> x;
    std::size_t k{};
    UIntType carry;
};

using ranlux24_base =
    subtract_with_carry_engine<std::uint_fast32_t, 24, 10, 24>;

using ranlux48_base = subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12>;
using ranlux48_base_correct =
    subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12, true>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_MERSENNE_TWISTER_ENGINE
