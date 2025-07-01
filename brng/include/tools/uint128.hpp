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

#ifndef BRNG_TOOLS_UINT128
#define BRNG_TOOLS_UINT128

#include <cstdlib>

#define FORCE_ENABLE 0

#if (!defined(__GNUC__) && !defined(__clang__)) || FORCE_ENABLE

#include <cstddef>
#include <cstdint>

namespace adhoc {

struct uint128 {
    std::uint64_t high;
    std::uint64_t low;

    constexpr uint128() : high(0), low(0) {}

    template <class T>
    constexpr uint128(T const v)
        : high(0), low(static_cast<std::uint64_t>(v)) {}

    template <class T>
    constexpr uint128(T h, T l)
        : high(static_cast<std::uint64_t>(h)),
          low(static_cast<std::uint64_t>(l)) {}

    explicit constexpr operator bool() const noexcept {
        return static_cast<bool>(this->high | this->low);
    }

    template <class T> constexpr operator T() const noexcept {
        return static_cast<T>(this->low);
    }

    auto constexpr operator==(const uint128 &v) const noexcept -> bool {
        return (this->low == v.low) && (this->high == v.high);
    }

    template <class T>
    auto constexpr operator==(T amount) const noexcept -> bool {
        return static_cast<std::uint64_t>(amount) == this->low;
    }

    auto constexpr operator+=(const uint128 &v) noexcept -> uint128 {
        auto o = this->low;
        this->low += v.low;
        this->high += v.high;
        this->high += static_cast<std::uint64_t>(this->low < o);
        return *this;
    }

    auto constexpr operator+(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t += v;
        return t;
    }

    template <class T> auto constexpr operator+(T v) const noexcept -> uint128 {
        uint128 t(*this);
        t += uint128(v);
        return t;
    }

    auto constexpr operator++() noexcept -> uint128 & {
        ++this->low;
        if (this->low == 0) {
            ++this->high;
        }
        return *this;
    }

    auto constexpr operator-=(const uint128 &v) noexcept -> uint128 {
        auto o = this->low;
        this->low -= v.low;
        this->high -= v.high;
        this->high -= static_cast<std::uint64_t>(this->low > o);
        return *this;
    }

    auto constexpr operator-(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t -= v;
        return t;
    }

    template <class T> auto constexpr operator-(T v) const noexcept -> uint128 {
        uint128 t(*this);
        t -= uint128(v);
        return t;
    }

    auto constexpr operator--() noexcept -> uint128 & {
        if (this->low == 0) {
            --this->high;
        }
        --this->low;
        return *this;
    }

    template <class T>
    auto constexpr operator<<(T amount) const noexcept -> uint128 {
        // uint64_t shifts of >= 64 are undefined, so we will need some
        // special-casing.
        return amount >= 64  ? uint128(this->low << (amount - 64U),
                                       static_cast<std::uint64_t>(0))
               : amount == 0 ? *this
                             : uint128((this->high << amount) |
                                           (this->low >> (64U - amount)),
                                       this->low << amount);
    }

    auto constexpr operator<<=(uint8_t amount) noexcept -> uint128 {
        if (amount == 0) {
            return *this;
        } else if (amount >= 128) {
            this->high = 0;
            this->low = 0;
            return *this;
        } else if (amount >= 64) {
            this->high = this->low << (amount - 64U);
            this->low = 0;
            return *this;
        }
        this->high <<= amount;
        this->high |= this->low >> (64U - amount);
        this->low <<= amount;
        return *this;
    }

    auto constexpr operator*(const uint128 &rhs) const noexcept -> uint128 {
        uint64_t a32 = this->low >> 32U;
        uint64_t a00 = this->low & 0xffffffff;
        uint64_t b32 = rhs.low >> 32U;
        uint64_t b00 = rhs.low & 0xffffffff;
        uint128 result = uint128(
            this->high * rhs.low + this->low * rhs.high + a32 * b32, a00 * b00);
        result += uint128(a32 * b00) << 32U;
        result += uint128(a00 * b32) << 32U;
        return result;
    }

    auto constexpr operator*=(const uint128 &rhs) noexcept -> uint128 {
        uint64_t a32 = this->low >> 32U;
        uint64_t a00 = this->low & 0xffffffff;
        uint64_t b32 = rhs.low >> 32U;
        uint64_t b00 = rhs.low & 0xffffffff;
        this->high = this->high * rhs.low + this->low * rhs.high + a32 * b32;
        this->low = a00 * b00;
        this->operator+=(uint128(a32 * b00) << 32U);
        this->operator+=(uint128(a00 * b32) << 32U);
        return *this;
    }

    auto constexpr operator<=(const uint128 &o) const noexcept -> bool {
        return this->high < o.high ||
               (this->high == o.high && this->low <= o.low);
    }

    auto constexpr operator>=(const uint128 &o) const noexcept -> bool {
        return this->high > o.high ||
               (this->high == o.high && this->low >= o.low);
    }

    template <class T>
    auto constexpr operator>>=(T amount) noexcept -> uint128 {
        static_assert(std::is_integral_v<T>);

        if (amount >= 64) {
            this->low = this->high >> (amount - 64);
            this->high = 0;
        } else if (amount != 0) {
            this->low = (this->high << (64 - amount)) | (this->low >> amount);
            this->high >>= amount;
        }
        return *this;
    }

    template <class T>
    auto constexpr operator>>(T v) const noexcept -> uint128 {
        uint128 t(*this);
        t >>= v;
        return t;
    }

    auto constexpr operator%=(const uint128 &b) -> uint128 {
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

    auto constexpr operator%(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t %= v;
        return t;
    }

    auto constexpr operator&=(const uint128 &b) -> uint128 {
        this->low &= b.low;
        this->high &= b.high;
        return *this;
    }

    auto constexpr operator&(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t &= v;
        return t;
    }

    template <class T> auto constexpr operator&(T v) const noexcept -> uint128 {
        uint128 t(*this);
        t.low &= static_cast<std::uint64_t>(v);
        return t;
    }

    auto constexpr operator^=(const uint128 &b) -> uint128 {
        this->low ^= b.low;
        this->high ^= b.high;
        return *this;
    }

    auto constexpr operator^(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t ^= v;
        return t;
    }

    template <class T> auto constexpr operator^(T v) const noexcept -> uint128 {
        uint128 t(*this);
        t.low ^= static_cast<std::uint64_t>(v);
        return t;
    }

    auto constexpr operator|=(const uint128 &b) -> uint128 {
        this->low |= b.low;
        this->high |= b.high;
        return *this;
    }

    auto constexpr operator|(const uint128 &v) const noexcept -> uint128 {
        uint128 t(*this);
        t |= v;
        return t;
    }

    template <class T> auto constexpr operator|(T v) const noexcept -> uint128 {
        uint128 t(*this);
        t.low |= static_cast<std::uint64_t>(v);
        return t;
    }

    auto constexpr operator~() const noexcept -> uint128 {
        uint128 t(*this);
        t.low = ~t.low;
        t.high = ~t.high;
        return t;
    }

    auto constexpr operator/(const uint128 &b) const noexcept -> uint128 {
        // if (!b)
        //     throw std::domain_error("divide by zero");

        uint128 x(b);
        uint128 y(*this >> 1U);
        unsigned int shift = 0;
        while (x <= y) {
            x <<= 1U;
            shift++;
        }

        uint128 result(0);
        uint128 thiscopy(*this);
        while (thiscopy >= b) {
            if (thiscopy >= x) {
                thiscopy -= x;
                result += (uint128{1} << shift);
            }
            shift--;
            x >>= 1U;
        }
        return result;
    }

    auto constexpr operator<(uint128 const &rhs) const -> bool {
        if (this->high < rhs.high) {
            return true;
        }

        if (this->high == rhs.high && this->low < rhs.low) {
            return true;
        }

        return false;
    }
};

struct div_t {
    uint128 rem;
    uint128 quot;
};

inline auto lldiv(uint128 const &x, uint128 const &y) noexcept -> div_t {
    // if (!y)
    //     throw std::domain_error("divide by zero");

    uint128 divisor(y);
    uint128 den_div2(x >> 1U);
    unsigned int shift = 0;
    while (divisor <= den_div2) {
        divisor <<= 1U;
        shift++;
    }

    div_t result{x, uint128(0)};
    while (result.rem >= y) {
        if (result.rem >= divisor) {
            result.rem -= divisor;
            result.quot += (uint128{1} << shift);
        }
        shift--;
        divisor >>= 1U;
    }

    return result;
}

auto constexpr max(uint128 const &lhs, uint128 const &rhs) -> uint128 {
    if (lhs < rhs) {
        return rhs;
    }
    return lhs;
}

} // namespace adhoc

template <> struct std::numeric_limits<adhoc::uint128> {
    static constexpr bool is_specialized = true;
    static constexpr int digits = 128;
    static constexpr adhoc::uint128 min() noexcept {
        return adhoc::uint128{0, 0};
    }
    static constexpr adhoc::uint128 max() noexcept {
        return adhoc::uint128{~0ULL, ~0ULL};
    }
    static constexpr bool is_signed = false;
    static constexpr bool is_integer = true;
    static constexpr bool is_exact = true;
    static constexpr int radix = 2;
};

namespace std {
template <> struct is_unsigned<adhoc::uint128> : std::true_type {};

// cland does not properly define is_unsigned_v. boo.
#if defined(__clang__)
template <>
constexpr bool is_unsigned_v<adhoc::uint128> =
    is_unsigned<adhoc::uint128>::value;
#endif

// MS compiler is even worse, it requires explicit true value for is_unsigned_v.
// double boo.
#if defined(_MSC_VER)
template <> constexpr bool is_unsigned_v<adhoc::uint128> = true;
#endif

} // namespace std

#else

namespace adhoc {

using uint128 = __uint128_t;

inline void lldiv() {}

} // namespace adhoc

#endif

constexpr auto operator""_ULLL(const char *digits) -> adhoc::uint128 {
    adhoc::uint128 result{};

    while (*digits != 0) {
        result *= 10;
        result += *digits - '0';
        ++digits;
    }

    return result;
}

#endif // BRNG_TOOLS_UINT128
