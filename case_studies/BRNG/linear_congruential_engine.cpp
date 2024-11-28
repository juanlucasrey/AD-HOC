#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "linear_congruential_engine.hpp"

// #include "double_integer.hpp"

#include <boost/multiprecision/cpp_int.hpp>
using boost::multiprecision::uint128_t;

#include <cstdint>
#include <limits>
#include <random>
#include <vector>

struct uint128 {
    std::uint64_t high;
    std::uint64_t low;

    template <class T>
    uint128(T v) : high(0), low(static_cast<std::uint64_t>(v)) {}

    template <class T>
    uint128(T h, T l)
        : high(static_cast<std::uint64_t>(h)),
          low(static_cast<std::uint64_t>(l)) {}

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

    auto operator*(uint128 rhs) noexcept -> uint128 {
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
};

int main() {
    TEST_START;

    // check against std minstd_rand0
    {
        std::minstd_rand0 rng1;
        adhoc::minstd_rand0 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // minstd_rand0 fwd and back
    {
        adhoc::minstd_rand0 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // minstd_rand0 back and fwd
    {
        adhoc::minstd_rand0 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std minstd_rand
    {
        std::minstd_rand rng1;
        adhoc::minstd_rand rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // minstd_rand fwd and back
    {
        adhoc::minstd_rand rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // minstd_rand back and fwd
    {
        adhoc::minstd_rand rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std ZX81
    {
        std::linear_congruential_engine<std::uint_fast16_t, 75, 74, 0> rng1;
        adhoc::ZX81 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // ZX81 fwd and back
    {
        adhoc::ZX81 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // ZX81 back and fwd
    {
        adhoc::ZX81 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std ranqd1
    {
        std::linear_congruential_engine<std::uint_fast32_t, 1664525, 1013904223,
                                        0>
            rng1;
        adhoc::ranqd1 rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // ranqd1 fwd and back
    {
        adhoc::ranqd1 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // ranqd1 back and fwd
    {
        adhoc::ranqd1 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std RANDU
    {
        std::linear_congruential_engine<std::uint_fast32_t, 65539, 0,
                                        2147483648>
            rng1;
        adhoc::RANDU rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // RANDU fwd and back
    {
        adhoc::RANDU rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // RANDU back and fwd
    {
        adhoc::RANDU rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // check against std Borland
    {
        std::linear_congruential_engine<std::uint_fast32_t, 22695477, 1,
                                        2147483648>
            rng1;
        adhoc::Borland rng2;
        TEST_FUNC(compare_rng(rng1, rng2, 100));
    }

    // Borland fwd and back
    {
        adhoc::Borland rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // Borland back and fwd
    {
        adhoc::Borland rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    std::uint64_t val1 =
        static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max());
    std::uint64_t val2 =
        static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max());
    std::uint64_t val3 = val1 * val2;
    std::cout << val3 << std::endl;

    {
        uint128 val11(std::numeric_limits<std::uint64_t>::max());
        uint128 val12(std::numeric_limits<std::uint64_t>::max());
        auto val13 = val11 * val12;
        std::cout << "stop" << std::endl;
    }

    // typedef double_integer<unsigned long long, unsigned long long> uint128_t;

    {
        uint128_t val11(std::numeric_limits<std::uint64_t>::max());
        uint128_t val12(std::numeric_limits<std::uint64_t>::max());
        uint128_t val13 = val11 * val12;
        std::cout << val13 << std::endl;
    }

    TEST_END;
}
