#include <tuple_utils.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

class type1 {};
class type2 {};
class type3 {};
class type4 {};
class type_free {};

TEST(TupleUtils, Remove) {
    constexpr auto tup =
        std::tuple<type4, type3, type2, type1, type_free, type_free, type_free,
                   type_free, type_free, type_free, type_free>{};

    constexpr auto bt_free = replace(tup, type4{}, type_free{});

    std::cout << type_name2<decltype(bt_free)>() << std::endl;
}

} // namespace adhoc3
