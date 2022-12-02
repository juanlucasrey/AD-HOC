#include <new/tuple_utils.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc2 {

class available_t {};
class replacement_t {};
class node1 {};
class node2 {};
class node3 {};
class node4 {};
class node5 {};
class node6 {};

TEST(TupleUitls, First) {
    auto constexpr res3 = replace_first<available_t, replacement_t>(
        std::tuple<available_t, node1, node2, node3, node4, node5, node6>{});

    static_assert(std::is_same_v<decltype(res3),
                                 const std::tuple<replacement_t, node1, node2,
                                                  node3, node4, node5, node6>>);
}

TEST(TupleUitls, Middle) {
    auto constexpr res3 = replace_first<available_t, replacement_t>(
        std::tuple<node1, node2, node3, available_t, node4, node5, node6>{});

    static_assert(
        std::is_same_v<decltype(res3),
                       const std::tuple<node1, node2, node3, replacement_t,
                                        node4, node5, node6>>);
}

TEST(TupleUitls, End) {
    auto constexpr res3 = replace_first<available_t, replacement_t>(
        std::tuple<node1, node2, node3, node4, node5, node6, available_t>{});

    static_assert(
        std::is_same_v<decltype(res3),
                       const std::tuple<node1, node2, node3, node4, node5,
                                        node6, replacement_t>>);
}

TEST(TupleUitls, Repeat) {
    auto constexpr res3 = replace_first<available_t, replacement_t>(
        std::tuple<available_t, available_t, node1, node2, node3, node4, node5,
                   node6>{});

    static_assert(
        std::is_same_v<decltype(res3),
                       const std::tuple<replacement_t, available_t, node1,
                                        node2, node3, node4, node5, node6>>);
}

TEST(TupleUitls, GenerateTuple) {
    auto constexpr temp = generate_tuple_type<available_t, 3>::type{};
    static_assert(std::is_same_v<
                  decltype(temp),
                  const std::tuple<available_t, available_t, available_t>>);
}

} // namespace adhoc2