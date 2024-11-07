#ifndef TEST_TOOLS_HPP
#define TEST_TOOLS_HPP

#include <cmath>
#include <iostream>
#include <limits>

template <class D>
inline std::tuple<bool, double> expect_near_abs(D val1, D val2, D tol) {
    const D abs_diff = std::abs(val1 - val2);
    return std::make_tuple(abs_diff < std::abs(tol), abs_diff);
}

template <class D>
inline std::tuple<bool, double> expect_near_rel(D val1, D val2, D tol) {

    // should take care of the val1=val2=0 case
    if (val1 == val2) {
        return std::make_tuple(true, 0.);
    }

    D average = (val1 + val2) * 0.5;

    if (average == 0.) {
        average = std::max(std::abs(val1), std::abs(val2));
    }

    D rel_diff = (val1 - val2) / average;

    return std::make_tuple(rel_diff < tol, rel_diff);
}

#define TEST_START int _result = 0
#define TEST_END return _result

#define EXPECT_NEAR_ABS(VAL1, VAL2, TOL)                                       \
    {                                                                          \
        auto [is_near, tol] = expect_near_abs(VAL1, VAL2, TOL);                \
        if (!is_near) {                                                        \
            std::cout.precision(std::numeric_limits<double>::max_digits10);    \
            std::cout << __FILE__ << ":" << __LINE__ << " Failure"             \
                      << std::endl;                                            \
            std::cout << "Absolute difference " << tol << " exceeds " << TOL   \
                      << ", where" << std::endl;                               \
            std::cout << "val1 evaluates to " << VAL1 << std::endl;            \
            std::cout << "val2 evaluates to " << VAL2 << std::endl;            \
            _result = 1;                                                       \
        }                                                                      \
    }

#define EXPECT_NEAR_REL(VAL1, VAL2, TOL)                                       \
    {                                                                          \
        auto [is_near, tol] = expect_near_rel(VAL1, VAL2, TOL);                \
        if (!is_near) {                                                        \
            std::cout.precision(std::numeric_limits<double>::max_digits10);    \
            std::cout << __FILE__ << ":" << __LINE__ << " Failure"             \
                      << std::endl;                                            \
            std::cout << "Relative difference " << tol << " exceeds " << TOL   \
                      << ", where" << std::endl;                               \
            std::cout << "val1 evaluates to " << VAL1 << std::endl;            \
            std::cout << "val2 evaluates to " << VAL2 << std::endl;            \
            _result = 1;                                                       \
        }                                                                      \
    }

#endif // TEST_TOOLS_HPP
