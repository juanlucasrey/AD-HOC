#include <partition/binomial_coefficient.hpp>

#include <gtest/gtest.h>

#include <utility>

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

TEST(Partition, BinomialCoefficients) {
    static_assert(std::is_same_v<decltype(BinomialCoefficients<0>()),
                                 std::index_sequence<1>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficients<1>()),
                                 std::index_sequence<1, 1>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficients<2>()),
                                 std::index_sequence<1, 2, 1>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficients<3>()),
                                 std::index_sequence<1, 3, 3, 1>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficients<4>()),
                                 std::index_sequence<1, 4, 6, 4, 1>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficients<5>()),
                                 std::index_sequence<1, 5, 10, 10, 5, 1>>);
}

TEST(Partition, BinomialCoefficientsSquare) {
    static_assert(std::is_same_v<decltype(BinomialCoefficientsSquare<0>()),
                                 std::index_sequence<1>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficientsSquare<1>()),
                                 std::index_sequence<2>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficientsSquare<2>()),
                                 std::index_sequence<2, 2>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficientsSquare<3>()),
                                 std::index_sequence<2, 6>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficientsSquare<4>()),
                                 std::index_sequence<2, 8, 6>>);

    static_assert(std::is_same_v<decltype(BinomialCoefficientsSquare<5>()),
                                 std::index_sequence<2, 10, 20>>);
}

} // namespace adhoc3
