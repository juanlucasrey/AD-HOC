#include <adhoc.hpp>
#include <differential_operator.hpp>
#include <name.hpp>

#include <gtest/gtest.h>

namespace adhoc4 {

namespace {

inline constexpr auto test1() {
    constexpr std::string_view tmp{"d(var)*d(var2*var)"};
    std::array<char, tmp.size()> res;
    std::copy_n(tmp.data(), tmp.size(), res.data());
    return res;
}

inline constexpr auto test2() {
    constexpr std::string_view tmp{"d^2(var)*d^2(var2*var)*d(cos(var))"};
    std::array<char, tmp.size()> res;
    std::copy_n(tmp.data(), tmp.size(), res.data());
    return res;
}

} // namespace

TEST(Name, Init) {
    ADHOC(var);
    ADHOC(var2);
    auto res = d(var) * d(var2 * var);

    constexpr auto varprint = name(res);
    static_assert(test1() == varprint);

    auto res2 = res * res * d(cos(var));
    constexpr auto varprint2 = name(res2);
    static_assert(test2() == varprint2);
}

} // namespace adhoc4
