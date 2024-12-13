#ifndef CASE_STUDIES_BRNG_SUBTRACT_WITH_CARRY_ENGINE
#define CASE_STUDIES_BRNG_SUBTRACT_WITH_CARRY_ENGINE

#include "linear_congruential_engine.hpp"

#include <algorithm>
#include <array>
#include <cmath>

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
          bool original = false>
class subtract_with_carry_engine final {
  private:
    static constexpr std::size_t k = std::size_t(w / 32) + 1;

    void init(std::array<std::uint_least32_t, r * k> const &seeds) {
        auto iter = seeds.begin();
        for (std::size_t j = 0; j < r; j++) {
            UIntType val = 0;
            for (std::size_t ki = 0; ki < k; ++ki) {
                val += static_cast<UIntType>(*iter++) << 32 * ki;
            }
            this->x[j] = val % modulus;
        }

        this->carry = (this->x[long_lag - 1] == 0);
        if constexpr (!original) {
            for (std::size_t j = 0; j < long_lag; ++j) {
                this->operator()();
            }
            for (std::size_t j = 0; j < long_lag; ++j) {
                this->operator()<false>();
            }
        }
    }

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t short_lag = s;
    static constexpr std::size_t long_lag = r;
    static constexpr std::uint_least32_t default_seed = 19780503U;

    static constexpr result_type modulus = static_cast<result_type>(1U) << w;

    static_assert(s >= 1U);
    static_assert(s < r);
    static_assert(w >= 1U);
    static_assert(w <= std::numeric_limits<UIntType>::digits);

    subtract_with_carry_engine() : subtract_with_carry_engine(0U) {}

    explicit subtract_with_carry_engine(result_type value) {
        linear_congruential_engine<std::uint_least32_t, 40014U, 0U, 2147483563U>
            e(value == 0U ? default_seed : value);
        std::array<std::uint_least32_t, r * k> seeds;
        std::generate(seeds.begin(), seeds.end(), e);
        this->init(seeds);
    }

    template <class SeedSeq> explicit subtract_with_carry_engine(SeedSeq &seq) {
        std::array<std::uint_least32_t, r * k> a;
        seq.generate(a.begin(), a.end());
        this->init(a);
    }

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        if constexpr (FwdDirection) {
            const std::size_t short_index =
                (this->i < short_lag) ? (this->i + long_lag - short_lag)
                                      : (this->i - short_lag);

            const UIntType carry_prev = carry;
            const UIntType x_prev = this->x[this->i];
            const UIntType temp = this->x[this->i] + carry;
            if (this->x[short_index] >= temp) {
                this->x[this->i] = this->x[short_index] - temp;
                this->carry = 0;
            } else {
                this->x[this->i] = modulus - temp + this->x[short_index];
                this->carry = 1;
            }

            const UIntType result = this->x[this->i];
            this->i = (this->i == (long_lag - 1)) ? 0 : (this->i + 1);
            return result;
        } else {
            this->i = (this->i == 0) ? (long_lag - 1) : (this->i - 1);
            const UIntType result = this->x[this->i];

            const std::size_t short_index =
                (this->i < short_lag) ? (this->i + long_lag - short_lag)
                                      : (this->i - short_lag);

            const UIntType temp =
                this->carry ? modulus - this->x[this->i] + this->x[short_index]
                            : this->x[short_index] - this->x[this->i];

            if (temp == 0) {
                this->carry = 0;
            } else if (temp == modulus) {
                this->carry = 1;
            } else {
                std::size_t k_prev = this->i;
                UIntType temp_prev = 0;
                std::size_t short_index_prev = short_index;
                do {
                    k_prev = (k_prev == 0) ? (long_lag - 1) : (k_prev - 1);
                    short_index_prev = (k_prev < short_lag)
                                           ? (k_prev + long_lag - short_lag)
                                           : (k_prev - short_lag);
                    temp_prev = this->x[k_prev];
                    if (temp_prev > this->x[short_index_prev]) {
                        temp_prev =
                            modulus - temp_prev + this->x[short_index_prev];
                    } else {
                        temp_prev = this->x[short_index_prev] - temp_prev;
                    }
                } while (temp_prev == 0 && k_prev != this->i);

                if (this->x[short_index_prev] >= temp_prev) {
                    this->carry = 0;
                } else {
                    this->carry = 1;
                }
            }

            this->x[this->i] = temp - this->carry;

            return result;
        }
    }

    void discard(unsigned long long z) {
        for (unsigned long long j = 0; j < z; ++j) {
            this->operator()();
        }
    }

    static constexpr auto min() -> UIntType {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> UIntType {
        constexpr UIntType ach = static_cast<UIntType>(1U)
                                 << static_cast<std::size_t>(w - 1U);

        return ach | (ach - static_cast<UIntType>(1U));

        // equivalent to the following but shift does not overflow
        // return (static_cast<result_type>(1U) << w) - 1U;
    }

    auto operator==(const subtract_with_carry_engine &rhs) const -> bool {
        return (this->carry == rhs.carry) && (this->i == rhs.i) &&
               (this->x == rhs.x);
    }

  private:
    std::array<UIntType, long_lag> x{0};
    std::size_t i{0};
    UIntType carry{0};
};

using ranlux24_base =
    subtract_with_carry_engine<std::uint_fast32_t, 24, 10, 24>;

using ranlux48_base = subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_MERSENNE_TWISTER_ENGINE
