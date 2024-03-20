#include "../../public3/partition/binomial_coefficient.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Partition, BinomialCoefficient) {
    // std::cout << "n=0" << std::endl;
    static_assert(binomial_coefficient(0, 0) == 1);
    // std::cout << "n=1" << std::endl;
    static_assert(binomial_coefficient(1, 0) == 1);
    static_assert(binomial_coefficient(1, 1) == 1);

    static_assert(binomial_coefficient(2, 0) == 1);
    static_assert(binomial_coefficient(2, 1) == 2);
    static_assert(binomial_coefficient(2, 2) == 1);
    static_assert(binomial_coefficient(3, 0) == 1);
    static_assert(binomial_coefficient(3, 1) == 3);
    static_assert(binomial_coefficient(3, 2) == 3);
    static_assert(binomial_coefficient(3, 3) == 1);
    // n=4
    static_assert(binomial_coefficient(4, 0) == 1);
    static_assert(binomial_coefficient(4, 1) == 4);
    static_assert(binomial_coefficient(4, 2) == 6);
    static_assert(binomial_coefficient(4, 3) == 4);
    static_assert(binomial_coefficient(4, 4) == 1);
}

} // namespace adhoc3
