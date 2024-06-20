#include <adhoc.hpp>
#include <base.hpp>
#include <name.hpp>

#include <gtest/gtest.h>

#include "type_name.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string_view>

namespace adhoc4 {

TEST(Adhoc, Types) {
    ADHOC(vol);
    ADHOC(strike);
    ADHOC(spot);

    auto res = vol * strike + exp(spot);

    auto res2 = comp_ellint_1(asin(exp(spot)));

    constexpr auto namestr = name(strike);
    std::cout << std::string_view(namestr.data(), namestr.size()) << std::endl;

    constexpr auto namestr2 = name(res2);
    std::cout << std::string_view(namestr2.data(), namestr2.size())
              << std::endl;
}

} // namespace adhoc4
