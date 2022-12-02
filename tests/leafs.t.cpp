#include <new/init.hpp>
#include <new/leafs.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(Leafs, First) {

    auto [v0, v1, v2, v3] = Init<4>();
    auto res1 = v0 * v0 * v1 + v2 * v3;
    auto res2 = v0 + v1 * v2 + v3;

    constexpr auto res = Leafs(res1, res2);
    const auto temp = std::make_tuple(decltype(v0){}, decltype(v1){},
                                      decltype(v2){}, decltype(v3){});
    static_assert(std::is_same_v<decltype(res), decltype(temp)>);
}

} // namespace adhoc2
