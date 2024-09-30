#include <adhoc.hpp>
#include <monomial_included.hpp>
#include <utils/flatten.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc4 {

TEST(MonomialIncluded, FlattenTest) {
    constexpr auto in = std::make_tuple(
        std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                        std::index_sequence<0, 1, 0, 0, 0, 0>{}),
        std::make_tuple(std::index_sequence<0, 0, 1, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 1, 0, 0>{}),
        std::make_tuple(std::index_sequence<0, 0, 0, 0, 1, 0>{},
                        std::index_sequence<0, 0, 0, 0, 0, 1>{}));

    constexpr auto out = flatten(in);

    constexpr auto out2 = std::make_tuple(
        std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 1, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 1, 0>{}),
        std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 1, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 0, 1>{}),
        std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 1, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 1, 0>{}),
        std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 1, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 0, 1>{}),
        std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 1, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 1, 0>{}),
        std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 1, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 0, 1>{}),
        std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 1, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 1, 0>{}),
        std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                        std::index_sequence<0, 0, 0, 1, 0, 0>{},
                        std::index_sequence<0, 0, 0, 0, 0, 1>{}));

    static_assert(std::is_same_v<decltype(out), decltype(out2)>);
}

TEST(MonomialIncluded, ErrorFix) {
    ADHOC(v1);
    auto res = erfc(v1);
    auto out1 = d<1>(res) * d(v1);
    {
        auto res = d(v1);
        auto res_t = std::make_tuple(res);
        constexpr auto resbool = monomial_included(out1, res_t);
        // EXPECT_FALSE(resbool);
        static_assert(!resbool);
    }

    {
        auto res = d<2>(v1);
        auto res_t = std::make_tuple(res);
        constexpr auto resbool = monomial_included(out1, res_t);
        // EXPECT_TRUE(resbool);
        static_assert(resbool);
    }

    {
        auto res3 = d<3>(v1);
        auto res_t = std::make_tuple(res3);
        constexpr auto resbool = monomial_included(out1, res_t);
        // EXPECT_TRUE(resbool);
        static_assert(resbool);
    }
}

TEST(MonomialIncluded, SumOrder2) {
    ADHOC(x);
    ADHOC(y);
    constexpr auto resbool = monomial_included(
        d<2>(x + y), std::make_tuple(d(x), d(y), d(x) * d(y)));

    // EXPECT_TRUE(resbool);
    static_assert(resbool);
}

} // namespace adhoc4
