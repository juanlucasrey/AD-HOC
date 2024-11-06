#ifndef TEST_TOOLS_HPP
#define TEST_TOOLS_HPP

#include <cmath>
#include <iostream>

template <class D> inline bool expect_near(D val1, D val2, D tol) {
    return std::abs(val1 - val2) < tol;
}

#define EXPECT_NEAR(VAL1, VAL2, TOL)                                           \
    if (!expect_near(VAL1, VAL2, TOL)) {                                       \
        std::cout << __FILE__ << ":" << __LINE__ << " Failure" << std::endl;   \
        std::cout << "Difference exceeds " << TOL << ", where" << std::endl;   \
        std::cout << "val1 evaluates to " << VAL1 << std::endl;                \
        std::cout << "val2 evaluates to " << VAL2 << std::endl;                \
        return 1;                                                              \
    }

#endif // TEST_TOOLS_HPP
