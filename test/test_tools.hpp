#ifndef TEST_TOOLS_HPP
#define TEST_TOOLS_HPP

#include <cmath>

template <class D> inline bool expect_near(D val1, D val2, D tol) {
    return std::abs(val1 - val2) < tol;
}

#endif // TEST_TOOLS_HPP
