#ifndef TEST_TOOLS_RNG_HPP
#define TEST_TOOLS_RNG_HPP

#include "../../test/test_tools.hpp"

#include <random>

template <class RNG1, class RNG2>
auto compare_rng(RNG1 &rng1, RNG2 &rng2, std::size_t n) -> int {
    TEST_START;

    for (std::size_t i = 0; i < n; ++i) {
        auto val1 = rng1();
        auto val2 = rng2();
        if (val1 != val2) {
            std::cout << "diff" << val1 << ", " << val2 << ", " << i
                      << std::endl;
        }
        EXPECT_EQUAL(rng1(), rng2());
    }

    TEST_END;
}

template <class RNG1, class RNG2>
auto compare_rng_limits(RNG1 & /* rng1 */, RNG2 & /* rng2 */) -> int {
    TEST_START;

    EXPECT_EQUAL(RNG1::min(), RNG2::min());
    EXPECT_EQUAL(RNG1::max(), RNG2::max());

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
