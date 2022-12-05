#include <calc_order.hpp>
#include <init.hpp>
#include <tape_size_bwd.hpp>
#include <tape_size_fwd.hpp>

#include "call_price.hpp"
#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc {

TEST(TapeSizeFwd, Zero) {
    double_t<0> val1;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(val1))>;
    static_assert(size == 0);
}

TEST(TapeSizeFwd, tapefwdbi) {
    double_t<0> val1;
    double_t<1> val2;
    auto valprod = val1 * val2;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(valprod))>;
    static_assert(size == 1);
}

TEST(TapeSizeFwd, tapefwddiv) {
    double_t<0> val1;
    double_t<1> val2;
    auto result = val1 / val2;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(result))>;
    static_assert(size == 1);
}

TEST(TapeSizeFwd, tapefwd2complicated) {
    double_t<0> val1;
    double_t<1> val2;
    double_t<2> val3;
    auto result = val1 * val2 + exp(val3);
    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(result))>;
    static_assert(size == 3);
}

TEST(TapeSizeFwd, tapefwd2skip) {
    double_t<0> val1;
    double_t<1> val2;
    double_t<2> val3;
    auto temp = val1 * val2;
    auto temp2 = temp * val3;
    auto result = temp * temp2;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(result))>;
    static_assert(size == 3);
}

TEST(TapeSizeFwd, tapesizefwd2) {
    double_t<0> val1;
    double_t<1> val2;
    double_t<2> val3;
    auto temp = val1 * val2;
    auto temp2 = temp * val3;
    auto result = temp * temp2;

    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(result))>;
    static_assert(size == 3);
}

TEST(TapeSizeFwd, tapesizefwdExternalNode) {
    double_t<0> val1;
    double_t<1> val2;
    double_t<2> val3;
    double_t<2> val4;
    auto temp = val1 * val2;
    auto temp2 = temp * val3;
    auto result = temp * temp2;

    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(result))>;
    static_assert(size == 3);
}

TEST(TapeSizeFwd, BSSinglePrice) {
    auto [S, K, v, T] = Init<4>();
    auto result_adhoc = call_price(S, K, v, T);

    auto constexpr temp = tape_size_bwd<decltype(S), decltype(K), decltype(v),
                                        decltype(T), decltype(result_adhoc)>();

    static_assert(temp == 3);
}

TEST(TapeSizeBwd, CutLeafs) {
    auto [v0, v1, v2, v3] = Init<4>();
    auto result = (v0 * v1) * (v2 * v3);
    auto constexpr temp =
        tape_size_bwd<decltype(v0), decltype(v1), decltype(v2), decltype(v3),
                      decltype(result)>();
    auto constexpr temp2 =
        tape_size_bwd<decltype(v0), decltype(v1), decltype(result)>();

    static_assert(temp == 2);
    static_assert(temp2 == 1);
}

} // namespace adhoc
