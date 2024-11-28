#ifndef CASE_STUDIES_BRNG_LINEAR_CONGRUENTIAL_ENGINE
#define CASE_STUDIES_BRNG_LINEAR_CONGRUENTIAL_ENGINE

#include <cstdint>
#include <limits>

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

template <class UIntType>
constexpr auto modular_multiplicative_inverse_max_plus_one(UIntType b) {
    UIntType t1 = 0;
    bool t1_b = true; // true means positive, false negative
    UIntType t2 = 1;
    bool t2_b = true; // true means positive, false negative
    UIntType t3 = 0;
    bool t3_b = true; // true means positive, false negative

    // awkward first step to avoid overflow
    UIntType n_minus_one = std::numeric_limits<UIntType>::max();
    UIntType r = n_minus_one - (b - 1);
    UIntType q = 1 + r / b;
    r = r % b;

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

    UIntType n = b;
    b = r;
    t1 = t2;
    t1_b = t2_b;
    t2 = t3;
    t2_b = t3_b;

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
        return t1;
    } else {
        UIntType q = n_minus_one - (t1 - 1);
        return q;
    }
}

struct uint128 {
    std::uint64_t high;
    std::uint64_t low;

    template <class T>
    constexpr uint128(T const v)
        : high(0), low(static_cast<std::uint64_t>(v)) {}

    template <class T>
    uint128(T h, T l)
        : high(static_cast<std::uint64_t>(h)),
          low(static_cast<std::uint64_t>(l)) {}

    operator std::uint64_t() const { return low; }

    auto operator+=(const uint128 &v) noexcept -> uint128 {
        auto o = low;
        low += v.low;
        high += v.high;
        high += (low < o);
        return *this;
    }

    auto operator+(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t += v;
        return t;
    }

    uint128 operator-=(const uint128 &v) noexcept {
        auto o = low;
        low -= v.low;
        high -= v.high;
        high -= (low > o);
        return *this;
    }

    auto operator-(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t -= v;
        return t;
    }

    auto operator<<(int amount) -> uint128 {
        // uint64_t shifts of >= 64 are undefined, so we will need some
        // special-casing.
        return amount >= 64 ? uint128(low << (amount - 64),
                                      static_cast<std::uint64_t>(0))
               : amount == 0
                   ? *this
                   : uint128((high << amount) | (low >> (64 - amount)),
                             low << amount);
    }

    uint128 operator<<=(uint8_t amount) noexcept {
        high <<= amount;
        high |= low >> (64 - amount);
        low <<= amount;
        return *this;
    }

    auto operator*(const uint128 &rhs) const noexcept -> uint128 {
        uint64_t a32 = low >> 32;
        uint64_t a00 = low & 0xffffffff;
        uint64_t b32 = rhs.low >> 32;
        uint64_t b00 = rhs.low & 0xffffffff;
        uint128 result =
            uint128(high * rhs.low + low * rhs.high + a32 * b32, a00 * b00);
        result += uint128(a32 * b00) << 32;
        result += uint128(a00 * b32) << 32;
        return result;
    }

    bool operator<=(const uint128 &o) const noexcept {
        return high < o.high || (high == o.high && low <= o.low);
    }

    bool operator>=(const uint128 &o) const noexcept {
        return high > o.high || (high == o.high && low >= o.low);
    }

    uint128 operator>>=(uint8_t v) noexcept {
        low >>= v;
        low |= high << (64 - v);
        high >>= v;
        return *this;
    }

    uint128 operator>>(uint8_t v) const noexcept {
        uint128 t(*this);
        t >>= v;
        return t;
    }

    uint128 operator%=(const uint128 &b) {
        // if (!b)
        //     throw std::domain_error("divide by zero");

        uint128 x(b), y(*this >> static_cast<uint8_t>(1));
        while (x <= y) {
            x <<= 1;
        }
        while (*this >= b) {
            if (*this >= x)
                *this -= x;
            x >>= 1;
        }
        return *this;
    }

    uint128 operator%(const uint128 &v) const noexcept {
        uint128 t(*this);
        t %= v;
        return t;
    }
};

} // namespace detail

template <class UIntType, UIntType a, UIntType c, UIntType m>
class linear_congruential_engine final {
    static_assert(m == 0 || a < m);
    static_assert(m == 0 || c < m);
    static_assert(m == 0 ? detail::gcd(static_cast<UIntType>(2U), a) == 1U
                         : detail::gcd(m, a) == 1U);

  public:
    using result_type = UIntType;

    static constexpr UIntType multiplier = a;
    static constexpr UIntType multiplier_inverse =
        m == 0 ? detail::modular_multiplicative_inverse_max_plus_one(a)
               : detail::modular_multiplicative_inverse(m, a);
    static constexpr UIntType increment = c;
    static constexpr UIntType modulus = m;
    static constexpr UIntType default_seed = 1U;

    explicit linear_congruential_engine(result_type value = default_seed)
        : x(value) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        if constexpr (m == 0) {
            if constexpr (FwdDirection) {
                this->x = a * this->x + c;
                return this->x;
            } else {
                const auto result = this->x;
                this->x = multiplier_inverse * (this->x - c);
                return result;
            }
        } else {
            constexpr auto digits = std::numeric_limits<UIntType>::digits;

            using upgraded_type = std::conditional_t<
                digits == 8, std::uint_fast16_t,
                std::conditional_t<
                    digits == 16, std::uint_fast32_t,
                    std::conditional_t<digits == 32, std::uint_fast64_t,
                                       detail::uint128>>>;

            constexpr auto c_64 = static_cast<upgraded_type>(c);
            constexpr auto m_64 = static_cast<upgraded_type>(m);

            if constexpr (FwdDirection) {
                constexpr auto a_64 = static_cast<upgraded_type>(a);

                this->x = static_cast<UIntType>(
                    (a_64 * static_cast<upgraded_type>(this->x) + c_64) % m_64);
                return this->x;
            } else {
                constexpr auto a_inv_64 =
                    static_cast<upgraded_type>(multiplier_inverse);

                const auto result = this->x;
                this->x = static_cast<UIntType>(
                    (a_inv_64 * (static_cast<upgraded_type>(this->x) - c_64)) %
                    m_64);
                return result;
            }
        }
    }

    static constexpr auto min() -> UIntType { return c == 0U ? 1U : 0U; }
    static constexpr auto max() -> UIntType {
        return m == 0 ? std::numeric_limits<result_type>::max() : m - 1U;
    }

  private:
    UIntType x;
};

using minstd_rand0 =
    linear_congruential_engine<std::uint_fast32_t, 16807, 0, 2147483647>;

using minstd_rand =
    linear_congruential_engine<std::uint_fast32_t, 48271, 0, 2147483647>;

using ZX81 = linear_congruential_engine<std::uint_fast16_t, 75, 74, 0>;

using ranqd1 =
    linear_congruential_engine<std::uint_fast32_t, 1664525, 1013904223, 0>;

using RANDU =
    linear_congruential_engine<std::uint_fast32_t, 65539, 0, 2147483648>;

using Borland =
    linear_congruential_engine<std::uint_fast32_t, 22695477, 1, 2147483648>;

using Newlib =
    linear_congruential_engine<std::uint_fast64_t, 6364136223846793005, 1,
                               9223372036854775808U>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_LINEAR_CONGRUENTIAL_ENGINE
