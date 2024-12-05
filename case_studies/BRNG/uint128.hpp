#ifndef CASE_STUDIES_BRNG_UINT128
#define CASE_STUDIES_BRNG_UINT128

#include <cstddef>
#include <cstdint>

namespace adhoc {

class uint128 {
  private:
    std::uint64_t high;
    std::uint64_t low;

  public:
    template <class T>
    constexpr explicit uint128(T const v)
        : high(0), low(static_cast<std::uint64_t>(v)) {}

    template <class T>
    uint128(T h, T l)
        : high(static_cast<std::uint64_t>(h)),
          low(static_cast<std::uint64_t>(l)) {}

    explicit operator std::uint64_t() const { return low; }

    auto operator+=(const uint128 &v) noexcept -> uint128 {
        auto o = low;
        low += v.low;
        high += v.high;
        high += static_cast<std::uint64_t>(low < o);
        return *this;
    }

    auto operator+(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t += v;
        return t;
    }

    auto operator-=(const uint128 &v) noexcept -> uint128 {
        auto o = low;
        low -= v.low;
        high -= v.high;
        high -= static_cast<std::uint64_t>(low > o);
        return *this;
    }

    auto operator-(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t -= v;
        return t;
    }

    auto operator<<(unsigned int amount) const noexcept -> uint128 {
        // uint64_t shifts of >= 64 are undefined, so we will need some
        // special-casing.
        return amount >= 64 ? uint128(low << (amount - 64U),
                                      static_cast<std::uint64_t>(0))
               : amount == 0
                   ? *this
                   : uint128((high << amount) | (low >> (64U - amount)),
                             low << amount);
    }

    auto operator<<=(uint8_t amount) noexcept -> uint128 {
        high <<= amount;
        high |= low >> (64U - amount);
        low <<= amount;
        return *this;
    }

    auto operator*(const uint128 &rhs) const noexcept -> uint128 {
        uint64_t a32 = low >> 32U;
        uint64_t a00 = low & 0xffffffff;
        uint64_t b32 = rhs.low >> 32U;
        uint64_t b00 = rhs.low & 0xffffffff;
        uint128 result =
            uint128(high * rhs.low + low * rhs.high + a32 * b32, a00 * b00);
        result += uint128(a32 * b00) << 32U;
        result += uint128(a00 * b32) << 32U;
        return result;
    }

    auto operator<=(const uint128 &o) const noexcept -> bool {
        return high < o.high || (high == o.high && low <= o.low);
    }

    auto operator>=(const uint128 &o) const noexcept -> bool {
        return high > o.high || (high == o.high && low >= o.low);
    }

    auto operator>>(std::size_t v) const noexcept -> uint128 {
        uint128 t(*this);
        t >>= v;
        return t;
    }

    auto operator>>=(std::size_t amount) noexcept -> uint128 {
        if (amount >= 64) {
            low = high >> (amount - 64U);
            high = 0;
        } else {
            low = (high << (64U - amount)) | (low >> amount);
            high >>= amount;
        }
        return *this;
    }

    auto operator%=(const uint128 &b) -> uint128 {
        // if (!b)
        //     throw std::domain_error("divide by zero");

        uint128 x(b);
        uint128 y(*this >> 1U);
        while (x <= y) {
            x <<= 1U;
        }
        while (*this >= b) {
            if (*this >= x) {
                *this -= x;
            }
            x >>= 1U;
        }
        return *this;
    }

    auto operator%(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t %= v;
        return t;
    }
};

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_UINT128
