#ifndef TEST_TOOLS_RNG_HPP
#define TEST_TOOLS_RNG_HPP

#include "../../test/test_tools.hpp"

#include <random>

template <class RNG1, class RNG2>
auto compare_rng(RNG1 &rng1, RNG2 &rng2, std::size_t n) -> int {
    TEST_START;
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQUAL(rng1(), rng2());
    }

    std::uniform_real_distribution<> dist_r(1.0, 2.0);
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQUAL(dist_r(rng1), dist_r(rng2));
    }

    std::uniform_int_distribution<> dist_i(1, 6);
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQUAL(dist_i(rng1), dist_i(rng2));
    }
    TEST_END;
}

template <class RNG> auto check_fwd_and_back(RNG &rng, std::size_t n) -> int {
    TEST_START;
    std::vector<typename RNG::result_type> vals1(n);

    for (std::size_t i = 0; i < n; ++i) {
        vals1[i] = rng();
    }

    std::vector<typename RNG::result_type> vals2(n);
    for (std::size_t i = n; i--;) {
        vals2[i] = rng.template operator()<false>();
    }

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQUAL(vals1[i], vals2[i]);
    }
    TEST_END;
}

template <class RNG> auto check_back_and_fwd(RNG &rng, std::size_t n) -> int {
    TEST_START;
    std::vector<typename RNG::result_type> vals1(n);

    for (std::size_t i = 0; i < n; ++i) {
        vals1[i] = rng.template operator()<false>();
    }

    std::vector<typename RNG::result_type> vals2(n);
    for (std::size_t i = n; i--;) {
        vals2[i] = rng();
    }

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQUAL(vals1[i], vals2[i]);
    }
    TEST_END;
}

#endif // TEST_TOOLS_RNG_HPP
