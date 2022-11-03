#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>

#include <constants_double.hpp>

#include <gtest/gtest.h>

template <class T> constexpr std::string_view type_name() {
    using namespace std;
#ifdef __clang__
    string_view p = __PRETTY_FUNCTION__;
    return string_view(p.data() + 34, p.size() - 34 - 1);
#elif defined(__GNUC__)
    string_view p = __PRETTY_FUNCTION__;
#if __cplusplus < 201402
    return string_view(p.data() + 36, p.size() - 36 - 1);
#else
    return string_view(p.data() + 49, p.find(';', 49) - 49);
#endif
#elif defined(_MSC_VER)
    string_view p = __FUNCSIG__;
    return string_view(p.data() + 84, p.size() - 84 - 7);
#endif
}

namespace detail {
template <typename T> struct AsTemplateArg {
    T buffer{};
    constexpr AsTemplateArg(T t) : buffer(t) {}
    constexpr operator T() const { return this->buffer; }
};
} // namespace detail

template <typename T, detail::AsTemplateArg<T> F> struct C {
    constexpr operator T() const { return F; }
};

template <detail::AsTemplateArg<float> F> using CF = C<float, F>;
template <detail::AsTemplateArg<double> F> using CD = C<double, F>;
template <detail::AsTemplateArg<long double> F> using CLD = C<long double, F>;

TEST(adhoc2, nttp_final) {
    float constexpr onef = 1.0f;
    C<float, 1.0f> onef_type;
    // different sizes, type is empty
    static_assert(sizeof(onef) == 4);
    static_assert(sizeof(onef_type) == 1);
    // same values
    static_assert(onef == onef_type);

    double constexpr oned = 1.0;
    C<double, 1.0> oned_type;
    // different sizes, type is empty
    static_assert(sizeof(oned) == 8);
    static_assert(sizeof(oned_type) == 1);
    // same values
    static_assert(oned == oned_type);

    long double constexpr oneld = 1.0l;
    C<long double, 1.0l> oneld_type;
    // different sizes, type is empty
    static_assert(sizeof(oneld) == 16);
    static_assert(sizeof(oneld_type) == 1);
    // same values
    static_assert(oneld == oneld_type);

    double constexpr one = 1.0;
    CD<one> one_type;
    CD<one> one_type2;
    // same type for same values
    static_assert(std::is_same_v<decltype(one_type), decltype(one_type)>);
    static_assert(one == one_type);
    static_assert(one == one_type2);

    // 1.000000000000000222 = std::nextafter(1)
    double constexpr oneeps = 1.000000000000000222;
    CD<oneeps> oneeps_type;
    // different types for epsilon
    static_assert(!std::is_same_v<decltype(one_type), decltype(oneeps_type)>);
    static_assert(oneeps == oneeps_type);
}
