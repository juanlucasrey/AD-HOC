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

template <typename T> struct AsTemplateArg {
    T buffer{};
    constexpr AsTemplateArg(T t) : buffer(t) {}
    constexpr operator T() const { return this->buffer; }
};

template <typename T, AsTemplateArg<T> F> struct C {
    constexpr operator T() const { return F; }
};

TEST(adhoc2, nttp_final) {
    C<double, 0.4> tests7;
    C<float, 0.6f> tests8;
    C<long double, 0.4l> tests9;

    double constexpr one = 1.0;
    C<double, one> one_type;
    C<double, one> one_type2;
    static_assert(std::is_same_v<decltype(one_type), decltype(one_type)>);
    static_assert(one == one_type);
    double constexpr oneeps = 1.000000000000000222;
    C<double, oneeps> oneeps_type;
    static_assert(!std::is_same_v<decltype(one_type), decltype(oneeps_type)>);
    static_assert(oneeps == oneeps_type);

    std::cout << sizeof(one) << std::endl;
    std::cout << sizeof(one_type) << std::endl;
}
