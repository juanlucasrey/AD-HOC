#include <tape_size_fwd.hpp>

#include <gtest/gtest.h>

namespace adhoc {

TEST(TapeSizeFwd, Zero) {
    detail::adouble_aux<0> val1;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(fwd_calc_order_t(val1))>;
    static_assert(size == 0);
}

TEST(adhoc2, tapefwdbi) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    auto valprod = val1 * val2;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(fwd_calc_order_t(valprod))>;
    static_assert(size == 1);
}

TEST(adhoc2, tapefwddiv) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    auto result = val1 / val2;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(fwd_calc_order_t(result))>;
    static_assert(size == 1);
}

TEST(adhoc2, tapefwd2complicated) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    detail::adouble_aux<2> val3;
    auto result = val1 * val2 + exp(val3);
    constexpr std::size_t size =
        std::tuple_size_v<decltype(fwd_calc_order_t(result))>;
    static_assert(size == 3);
}

TEST(adhoc2, tapefwd2skip) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    detail::adouble_aux<2> val3;
    auto temp = val1 * val2;
    auto temp2 = temp * val3;
    auto result = temp * temp2;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(fwd_calc_order_t(result))>;
    static_assert(size == 3);
}

} // namespace adhoc
