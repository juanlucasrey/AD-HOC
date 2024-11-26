#ifndef CASE_STUDIES_BRNG_LINEAR_CONGRUENTIAL_ENGINE
#define CASE_STUDIES_BRNG_LINEAR_CONGRUENTIAL_ENGINE

#include <cstdint>

namespace adhoc {

namespace detail {

template <class UIntType>
constexpr auto gcd(UIntType a, UIntType b) -> UIntType {
    while (b != 0) {
        auto tmp = b;
        b = a % b;
        a = tmp;
    }
    return a;
};

template <class UIntType>
constexpr auto modular_multiplicative_inverse(UIntType n, UIntType b)
    -> UIntType {
    UIntType n_original = n;
    UIntType t1 = 0;
    bool t1_b = true; // true means positive, false negative
    UIntType t2 = 1;
    bool t2_b = true; // true means positive, false negative
    UIntType t3 = 0;
    bool t3_b = true; // true means positive, false negative
    UIntType r = 1;
    while (r != 0) {
        UIntType q = n / b;
        r = n % b;

        // t3 = t1 - q × t2
        if (t1_b && !t2_b) {
            t3_b = true;
            t3 = t1 + q * t2;
        } else if (!t1_b && t2_b) {
            t3_b = false;
            t3 = t1 + q * t2;
        } else if (t1_b && t2_b) {
            t3_b = t1 >= q * t2;
            if (t3_b) {
                t3 = t1 - q * t2;
            } else {
                t3 = q * t2 - t1;
            }
        } else {
            t3_b = q * t2 >= t1;
            if (t3_b) {
                t3 = q * t2 - t1;
            } else {
                t3 = t1 - q * t2;
            }
        }

        n = b;
        b = r;
        t1 = t2;
        t1_b = t2_b;
        t2 = t3;
        t2_b = t3_b;
    }

    if (t1_b) {
        UIntType q = t1 / n_original;
        r = t1 % n_original;
        return r;
    } else {
        UIntType q = t1 / n_original;
        r = t1 % n_original;
        return n_original - r;
    }
}

} // namespace detail

template <class UIntType, UIntType a, UIntType c, UIntType m>
class linear_congruential_engine final {
    static_assert(m == 0 || a < m);
    static_assert(m == 0 || c < m);
    static_assert(detail::gcd(m, a) == 1U);

  public:
    using result_type = UIntType;

    static constexpr UIntType multiplier = a;
    static constexpr UIntType increment = c;
    static constexpr UIntType modulus = m;
    static constexpr UIntType default_seed = 1U;

    explicit linear_congruential_engine(result_type value = default_seed)
        : x(value) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        constexpr auto a_64 = static_cast<std::uint_fast64_t>(a);
        constexpr auto c_64 = static_cast<std::uint_fast64_t>(c);
        constexpr auto m_64 = static_cast<std::uint_fast64_t>(m);

        if constexpr (FwdDirection) {
            this->x = static_cast<UIntType>(
                (a_64 * static_cast<std::uint_fast64_t>(this->x) + c_64) %
                m_64);
            return this->x;
        } else {
            constexpr auto a_inv = static_cast<std::uint_fast64_t>(
                detail::modular_multiplicative_inverse(m, a));

            const auto result = this->x;
            this->x = static_cast<UIntType>(
                (a_inv * (static_cast<std::uint_fast64_t>(this->x) - c_64)) %
                m_64);
            return result;
        }
    }

    static constexpr auto min() -> UIntType { return c == 0U ? 1U : 0U; }
    static constexpr auto max() -> UIntType { return m - 1U; }

  private:
    UIntType x;
};

using minstd_rand0 =
    linear_congruential_engine<std::uint_fast32_t, 16807, 0, 2147483647>;

using minstd_rand =
    linear_congruential_engine<std::uint_fast32_t, 48271, 0, 2147483647>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_LINEAR_CONGRUENTIAL_ENGINE
