#include "type_name.hpp"
#include <differential_operator.hpp>
#include <init.hpp>
#include <tape.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

template <class> struct is_tuple : std::false_type {};
template <class... T> struct is_tuple<std::tuple<T...>> : std::true_type {};

template <class... OutputsAndDerivatives>
constexpr auto Output2(OutputsAndDerivatives... /* in */) {
    return std::tuple<OutputsAndDerivatives...>{};
};

TEST(TapeInit, First) {
    auto [x, y] = Init<2>();
    auto var1 = x * y;
    auto var2 = cos(x * y);

    auto temp = Output2(var1, var2, d(x), d<2>(x) * d<2>(y));

    std::cout << "original tuple" << std::endl;
    std::cout << type_name2<decltype(temp)>() << std::endl;

    constexpr auto derivatives = std::apply(
        [](auto... ts) {
            return std::tuple_cat(
                std::conditional_t<is_tuple<decltype(ts)>::value,
                                   std::tuple<decltype(ts)>,
                                   std::tuple<>>{}...);
        },
        temp);

    std::cout << "derivatives" << std::endl;
    std::cout << type_name2<decltype(derivatives)>() << std::endl;

    constexpr auto variables = std::apply(
        [](auto... ts) {
            return std::tuple_cat(
                std::conditional_t<is_tuple<decltype(ts)>::value, std::tuple<>,
                                   std::tuple<decltype(ts)>>{}...);
        },
        temp);

    std::cout << "variables" << std::endl;
    std::cout << type_name2<decltype(variables)>() << std::endl;

    // Tape t(r, d(v0), d(v1));
    // t.set(v0) = 2.0;
    // t.set(v1) = 3.0;
    // t.der(r) = 1.0;

    // t.evaluate_fwd();
    // t.evaluate_bwd();

    // EXPECT_NEAR(t.der(v0), 1.0, 1e-13);
    // EXPECT_NEAR(t.der(v1), 1.0, 1e-13);
}

} // namespace adhoc3
