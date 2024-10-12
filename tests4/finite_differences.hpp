#ifndef TESTS4_FINITE_DIFFERENCES_HPP
#define TESTS4_FINITE_DIFFERENCES_HPP

#include <array>
#include <functional>

namespace adhoc4 {

std::array<double, 6> finite_differences(double x, double epsilon,
                                         std::function<double(double)> func);

} // namespace adhoc4

#endif // TESTS4_FINITE_DIFFERENCES_HPP
