#include <tape_size_fwd.hpp>

#include <gtest/gtest.h>

namespace adhoc {

TEST(adhoc2, tapefwd2complicated) {
    detail::double_t<0> val1;
    detail::double_t<1> val2;
    detail::double_t<2> val3;
    auto result = val1 * val2 + exp(val3);
    constexpr std::size_t size =
        std::tuple_size_v<decltype(fwd_calc_order_t(result))>;
    static_assert(size == 3);
}

TEST(adhoc2, tapefwd2skip) {
    detail::double_t<0> val1;
    detail::double_t<1> val2;
    detail::double_t<2> val3;
    auto temp = val1 * val2;
    auto temp2 = temp * val3;
    auto result = temp * temp2;
    constexpr std::size_t size =
        std::tuple_size_v<decltype(fwd_calc_order_t(result))>;
    static_assert(size == 3);
}

} // namespace adhoc
