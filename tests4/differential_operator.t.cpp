#include <adhoc.hpp>
#include <differential_operator.hpp>

#include <gtest/gtest.h>

namespace adhoc4 {

TEST(DifferentialOperator, Init) {
    ADHOC(var);
    ADHOC(var2);
    auto res = d(var) * d(var2 * var);

    static_assert(std::same_as<decltype(res * res),
                               decltype(d<2>(var) * d<2>(var2 * var))>);
}

} // namespace adhoc4
