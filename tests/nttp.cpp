#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>

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

template <typename T> void Ref(T &&);

typedef std::array<char, sizeof(float)> float_buf;

//////////////// ignore this boilerplate ////////////////
namespace std {

// In June 16, 2019 gcc, header <bit> exists but does not have bit_cast in it.
// stolen from https://en.cppreference.com/w/cpp/numeric/bit_cast
template <class To, class From>
typename std::enable_if<
    (sizeof(To) == sizeof(From)) && std::is_trivially_copyable<From>::value &&
        std::is_trivial<To>::value && true,
    // this implementation requires that To is trivially default constructible
    To>::type constexpr bit_cast(const From &src) noexcept {
    To dst = {};
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

// In June 16, 2019 gcc, bit_cast above is not really constexpr, so we have to
// define a custom bit_cast that works for 32-bit float *and* is constexpr.
template <>
constexpr float_buf bit_cast<float_buf, float>(const float &src) noexcept {
    if (src == 0)
        return {0, 0, 0, 0};

    bool sign = src < 0.0f;
    float f = sign ? -src : src;
    uint8_t exponent = 127;
    while (f >= 2 && exponent < 255) {
        exponent += 1, f *= 0.5f;
    }
    while (f < 1 && exponent > 1)
        exponent -= 1, f *= 2.0f;
    if (f < 1)
        exponent -= 1; // denormal
    uint32_t mantissa =
        exponent == 255 ? 0 : 0x00800000 * static_cast<std::uint32_t>(f);
    uint32_t rep = (sign ? 0x80000000 : 0U) +
                   static_cast<std::uint32_t>(exponent << 23U) +
                   static_cast<std::uint32_t>(mantissa & 0x7FFFFFU);
    return {(char)(rep), (char)(rep >> 8), (char)(rep >> 16),
            (char)(rep >> 24)};
}

} // namespace std
//////////////// end of boilerplate ////////////////

template <typename T> struct AsTemplateArg {
    std::array<char, sizeof(T)> buffer = {};
    constexpr AsTemplateArg(const std::array<char, sizeof(T)> buf)
        : buffer(buf) {}
    constexpr AsTemplateArg(T t)
        : AsTemplateArg(std::bit_cast<std::array<char, sizeof(T)>>(t)) {}
    constexpr operator T() const { return std::bit_cast<T>(this->buffer); }
};

template <AsTemplateArg<float> F> float float_doubler() { return float{F} * 2; }

template <AsTemplateArg<float> exponent> constexpr float pow(float base) {
    return exp(log(base) * float{exponent});
}

template <> constexpr float pow<AsTemplateArg<float>{1.0}>(float base) {
    return base;
}

template <> constexpr float pow<AsTemplateArg<float>{2.0}>(float base) {
    return base * base;
}

template <AsTemplateArg<float> exponent> constexpr auto type(float base) {
    return AsTemplateArg<float>{base};
}

template <class In1, AsTemplateArg<float> F> class Sum {};

template <class Derived> struct Base {
    using is_adhoc_tag = void;

    template <AsTemplateArg<float> exponent>
    constexpr auto operator+(float in) const -> auto;
};

template <class Derived>
template <AsTemplateArg<float> exponent>
constexpr auto Base<Derived>::operator+(float /* in */) const -> auto {
    return Sum<Derived, AsTemplateArg<float>{1.0}>();
}

class TreeCalc : public Base<TreeCalc> {};

constexpr auto double_to_uint64(double half) -> std::uint64_t {
    if (half > 0.5) {
        return 1;
    } else {
        return 2;
    }
}

struct Double {
    std::array<char, sizeof(double)> buffer = {};
    constexpr Double(const std::array<char, sizeof(double)> buf)
        : buffer(buf) {}
    constexpr Double(double t)
        : Double(std::bit_cast<std::array<char, sizeof(double)>>(t)) {}
    constexpr operator double() const {
        return std::bit_cast<double>(this->buffer);
    }
};

// template <std::uint64_t D> struct Double {
//     // std::uint64_t buffer = {};
//     std::array<char, sizeof(float)> buffer = {};
//     constexpr Double(const std::array<char, sizeof(float)> buf) : buffer(buf)
//     {} constexpr Double(float t)
//         : AsTemplateArg(std::bit_cast<std::array<char, sizeof(float)>>(t)) {}
// };

template <Double exponent> constexpr auto C(double base) {
    return Double{base};
}

template <typename T> struct AsTemplateArg3 {
    constexpr AsTemplateArg3(double) {}
    constexpr operator double() const { return 0.; }
};

template <AsTemplateArg<double> F, std::uint64_t Rep = double_to_uint64(F)>
struct C2 {};

TEST(adhoc2, nttp) {
    // #define demo(f)                                                                \
//     std::cout << "Twice " << (f) << " is " << float_doubler<f>() << '\n'
    //     demo(1.0);
    //     demo(2.37);
    //     demo(0.0);
    //     demo(-0.0);
    //     demo(INFINITY);
    //     demo(-INFINITY);
    //     demo(1e-37);
    //     demo(1.17549435E-38);
    //     demo(1.1754942E-38);
    //     demo(1.17E-38);
    //     demo(1e-39);
    //     demo(1e-44);

#define demo2(f) std::cout << (f) << " squared is " << pow<2.0f>(f) << '\n'
    demo2(7);
    demo2(0.0); // can't take the log of zero... but our pow works fine!
    std::cout << pow<2.0f>(7.0) << std::endl;
    std::cout << pow<1.0f>(7.0) << std::endl;
    std::cout << pow<2.1f>(7.0) << std::endl;

    auto constexpr temptype = type<2.0f>(3.0f);
    auto constexpr temp1 = pow<2.0f>(7.0);
    auto constexpr temp2 = pow<1.0f>(7.0);
    // auto constexpr temp3 = pow<1.1f>(7.0);

    std::cout << type_name<decltype(temptype)>() << std::endl;

    // TreeCalc some_calc;
    auto temp = C<2.0>;
    std::cout << type_name<decltype(temp)>() << std::endl;
    std::cout << sizeof(temp) << std::endl;
    // auto tempb = temp.buffer;
    // std::cout << temp.buffer[0] << std::endl;

    C2<0.4> tests;
    C2<0.6> tests2;
    std::cout << type_name<decltype(tests)>() << std::endl;
    std::cout << type_name<decltype(tests2)>() << std::endl;
    std::cout << sizeof(tests) << std::endl;
    std::cout << sizeof(tests2) << std::endl;
}
