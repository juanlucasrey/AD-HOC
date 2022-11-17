#include <new/calc_order.hpp>
#include <new/evaluate.hpp>
#include <new/init.hpp>
#include <new/tape_nodes.hpp>
#include <new/tape_size_fwd.hpp>

#include "../type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc2 {

// TEST(TapeSizeFwd, Zero) {
//     double_t<0> val1;
//     constexpr std::size_t size =
//         std::tuple_size_v<decltype(calc_order_t(val1))>;
//     static_assert(size == 0);
// }

// TEST(TapeSizeFwd, tapefwdbi) {
//     double_t<0> val1;
//     double_t<1> val2;
//     auto valprod = val1 * val2;
//     constexpr std::size_t size =
//         std::tuple_size_v<decltype(calc_order_t(valprod))>;
//     static_assert(size == 1);
// }

// TEST(TapeSizeFwd, tapefwddiv) {
//     double_t<0> val1;
//     double_t<1> val2;
//     auto result = val1 / val2;
//     constexpr std::size_t size =
//         std::tuple_size_v<decltype(calc_order_t(result))>;
//     static_assert(size == 1);
// }

// TEST(TapeSizeFwd, tapefwd2complicated) {
//     double_t<0> val1;
//     double_t<1> val2;
//     double_t<2> val3;
//     auto result = val1 * val2 + exp(val3);
//     constexpr std::size_t size =
//         std::tuple_size_v<decltype(calc_order_t(result))>;
//     static_assert(size == 3);
// }

// TEST(TapeSizeFwd, tapefwd2skip) {
//     double_t<0> val1;
//     double_t<1> val2;
//     double_t<2> val3;
//     auto temp = val1 * val2;
//     auto temp2 = temp * val3;
//     auto result = temp * temp2;
//     constexpr std::size_t size =
//         std::tuple_size_v<decltype(calc_order_t(result))>;
//     static_assert(size == 3);
// }

TEST(TapeSizeFwd, tapesizefwd) {
    double_t<0> val1;
    double_t<1> val2;
    double_t<2> val3;
    auto temp = val1 * val2;
    auto temp2 = temp * val3;
    auto result = temp * temp2;
    auto leaves_and_roots = TapeRootsAndLeafs(result);

    auto constexpr nodes_inside = tape_size_fwd_t(result);

    evaluate(leaves_and_roots);

    std::cout << type_name<decltype(nodes_inside)>() << std::endl;

    constexpr std::size_t size =
        std::tuple_size_v<decltype(calc_order_t(result))>;
    static_assert(size == 3);
}

template <class D> inline auto cdf_n(D const &x) {
    using constants::CD;
    using constants::encode;
    using std::erfc;
    constexpr double minus_one_over_root_two =
        -1.0 / constexpression::sqrt(2.0);
    return CD<encode(0.5)>() * erfc(x * CD<encode(minus_one_over_root_two)>());
    // return CD<0.5>() * erfc(x * CD<minus_one_over_root_two>());
}

template <class I1, class I2, class I3, class I4>
auto call_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using constants::encode;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<encode(0.5)>();
    auto d2 = d1 + totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

TEST(TapeSizeFwd, BSSinglePrice) {
    double result = 0.;

    {
        double S = 100.0;
        double K = 102.0;
        double v = 0.15;
        double T = 0.5;
        result = call_price(S, K, v, T);
    }

    auto [S, K, v, T] = Init<4>();
    auto result_adhoc = call_price(S, K, v, T);

    auto tape = TapeRootsAndLeafs(result_adhoc);
    tape.set(S, 100.0);
    tape.set(K, 102.0);
    tape.set(v, 0.15);
    tape.set(T, 0.5);

    evaluate(tape);
    double result2 = tape.get(result_adhoc);
    EXPECT_EQ(result2, result);

    auto int_tape = evaluate_fwd_return_vals(tape);
    auto totalvol = v * sqrt(T);
    double totalvol_val = int_tape.get(totalvol);
    auto totalvol_val2 = 0.15 * std::sqrt(0.5);
    EXPECT_EQ(totalvol_val, totalvol_val2);
}

} // namespace adhoc2
