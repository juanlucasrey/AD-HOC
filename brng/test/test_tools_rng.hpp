#ifndef BRNG_TEST_TOOLS_RNG
#define BRNG_TEST_TOOLS_RNG

#include "../../test_simple/test_simple_include.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

template <class RNG1, class RNG2>
void compare_rng(RNG1 &rng1, RNG2 &rng2, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        auto val1 = rng1();
        auto val2 = rng2();
        if (val1 != val2) {
            std::cout << "diff: " << static_cast<std::size_t>(val1) << ", "
                      << static_cast<std::size_t>(val2) << ", " << i
                      << std::endl;
        }
        EXPECT_EQUAL(val1, val2);
    }
}

template <class RNG1, class RNG2>
void compare_rng_limits(RNG1 & /* rng1 */, RNG2 & /* rng2 */) {
    EXPECT_EQUAL(RNG1::min(), RNG2::min());
    EXPECT_EQUAL(RNG1::max(), RNG2::max());
}

template <class RNG> void check_fwd_and_back(RNG &rng, std::size_t n) {
    std::vector<typename RNG::result_type> vals1(n);

    for (std::size_t i = 0; i < n - 1; ++i) {
        vals1[i] = rng();
    }
    vals1[n - 1] = rng();

    std::vector<typename RNG::result_type> vals2(n);
    vals2[n - 1] = rng.template operator()<false>();
    for (std::size_t i = n - 1; i--;) {
        vals2[i] = rng.template operator()<false>();
    }

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQUAL(vals1[i], vals2[i]);
    }
}

template <class RNG> void check_back_and_fwd(RNG &rng, std::size_t n) {
    std::vector<typename RNG::result_type> vals1(n);

    for (std::size_t i = 0; i < n - 1; ++i) {
        vals1[i] = rng.template operator()<false>();
    }
    vals1[n - 1] = rng.template operator()<false>();

    std::vector<typename RNG::result_type> vals2(n);
    vals2[n - 1] = rng();
    for (std::size_t i = n - 1; i--;) {
        vals2[i] = rng();
    }

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQUAL(vals1[i], vals2[i]);
    }
}

template <class RNG> void check_discard(RNG &rng, std::size_t n) {
    RNG rng2 = rng;
    for (std::size_t i = 0; i < n; ++i) {
        rng();
    }
    rng2.discard(n);
    EXPECT_EQUAL(rng, rng2);
}

#endif // BRNG_TEST_TOOLS_RNG
