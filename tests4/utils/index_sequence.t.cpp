#include <utils/index_sequence.hpp>

#include <gtest/gtest.h>

namespace adhoc4 {

TEST(IndexSequence, MakeIndexSequence) {
    {
        constexpr auto seq = make_index_sequence<3, 7>();
        constexpr auto res = std::index_sequence<3, 4, 5, 6, 7>{};

        static_assert(std::same_as<decltype(seq), decltype(res)>);
    }

    {
        constexpr auto seq = make_index_sequence<7, 3>();
        constexpr auto res = std::index_sequence<7, 6, 5, 4, 3>{};

        static_assert(std::same_as<decltype(seq), decltype(res)>);
    }
}

} // namespace adhoc4
