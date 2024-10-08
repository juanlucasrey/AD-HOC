
#include <utils/tuple.hpp>

#include "../type_name.hpp"
#include <gtest/gtest.h>

#include <tuple>

namespace adhoc4 {

TEST(UtilsTuple, ReplaceFirst) {

    constexpr std::tuple<double, int, char, int> temp;
    constexpr int v1 = 0;
    constexpr float v2 = 0;

    // constexpr auto first = head(temp);

    // constexpr auto newtuple = first_type_is2(temp, v1);
    constexpr auto newtuple = replace_first(temp, v1, v2);
    std::cout << type_name2<decltype(newtuple)>() << std::endl;

    constexpr auto newtuple2 = replace_first(temp, v2, v1);
    std::cout << type_name2<decltype(newtuple2)>() << std::endl;
}

TEST(UtilsTuple, Invert) {
    constexpr std::tuple<double, int, char> temp;
    std::cout << type_name2<decltype(temp)>() << std::endl;
    // double a = 0;
    // auto result = remove(temp, a);
    // std::cout << type_name2<decltype(result)>() << std::endl;

    constexpr auto inverted_temp = invert(temp);
    std::cout << type_name2<decltype(inverted_temp)>() << std::endl;
}

} // namespace adhoc4
