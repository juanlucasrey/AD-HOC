#include "../../../test_simple/test_simple_include.hpp"

#include <tools/mask.hpp>
#include <tools/modular_arithmetic.hpp>
#include <tools/uint128.hpp>

#include <mersenne_twister_engine.hpp>
#include <philox_engine.hpp>

#include <cstdint>

auto main() -> int {

    constexpr std::uint64_t a13n = 810728U;
    constexpr std::uint64_t a23n = 1370589U;
    constexpr std::uint64_t m1 = 4294967087U;
    constexpr std::uint64_t m2 = 4294944443U;

    constexpr auto a13n_inv = adhoc::modular_multiplicative_inverse(m1, a13n);
    static_assert(a13n_inv == 2349796154U);
    constexpr auto a23n_inv = adhoc::modular_multiplicative_inverse(m2, a23n);
    static_assert(a23n_inv == 69372715U);

#ifndef _MSC_VER
    {
        adhoc::philox_engine<adhoc::uint128, 128, 4, 10, 0xCA5A826395121157,
                             0x9E3779B97F4A7C15, 0xD2E7470EE14C6C93,
                             0xBB67AE8584CAA73B>
            rng;
        for (std::size_t i = 0; i < 10000; i++) {
            adhoc::uint128 val = rng();

            if (adhoc::gcd(val, adhoc::uint128(2)) != 1) {
                continue;
            }

            auto valinv = adhoc::modular_multiplicative_inverse_pow2(128, val);
            auto valinv2 =
                adhoc::modular_multiplicative_inverse_pow2<adhoc::uint128,
                                                           true>(128, val);
            EXPECT_EQUAL(valinv, valinv2);
        }
    }
#endif

    adhoc::mt19937_64 rng;

    for (std::size_t i = 0; i < 10000; i++) {
        std::uint64_t val = rng();
        adhoc::uint128 mod = 1;
        mod <<= 64;

        adhoc::uint128 val128 = val;

        if (adhoc::gcd(val128, mod) != 1) {
            continue;
        }

        auto valinv = adhoc::modular_multiplicative_inverse_pow2(64, val);
        auto valinv2 =
            adhoc::modular_multiplicative_inverse_pow2<std::uint64_t, true>(
                64, val);
        EXPECT_EQUAL(valinv, valinv2);

        adhoc::uint128 valinv128 = valinv;
        auto res = val128 * valinv128;

        res %= mod;

        EXPECT_EQUAL(res, 1);
    }

    for (std::size_t i = 0; i < 1000; i++) {
        std::uint64_t m1 = rng();
        for (std::size_t j = 0; j < 1000; j++) {

            std::uint64_t val = rng();
            val %= m1;

            if (adhoc::gcd(val, m1) != 1) {
                continue;
            }

            auto valinv = adhoc::modular_multiplicative_inverse(m1, val);

            adhoc::uint128 val128 = val;
            adhoc::uint128 valinv128 = valinv;
            auto res = val128 * valinv128;

            res %= m1;
            EXPECT_EQUAL(res, 1);
        }
    }

    TEST_END;
}
