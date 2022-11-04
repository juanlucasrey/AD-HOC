#include <adhoc2.hpp>
#include <tape_size.hpp>
#include <utils.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, tapesize) {
    adouble val1(1.);
    constexpr std::size_t size = tape_size(val1, val1);
    static_assert(size == 1);
}

TEST(adhoc2, tapesizezero) {
    adouble val1(1.);
    adouble val2(1.);
    constexpr std::size_t size = tape_size(val1, val2);
    static_assert(size == 0);
}

TEST(adhoc2, tapemult) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod = val1 * val2;
    constexpr std::size_t size = tape_size(valprod, val1, val2);
    static_assert(size == 2);
}

TEST(adhoc2, tapediv) {
    adouble val1(2.0);
    adouble val2(3.0);
    auto result = val1 / val2;
    constexpr std::size_t size = tape_size(result, val2);
    static_assert(size == 1);

    constexpr std::size_t size2 = tape_size(result, val1, val2);
    static_assert(size2 == 2);
}

TEST(adhoc2, tapeadd) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod = val1 + val2;
    constexpr std::size_t size = tape_size(valprod, val1, val2);
    static_assert(size == 2);
}

TEST(adhoc2, tapecomplexd) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2 + val3;
    auto valprod = val1 * val2;
    auto valprod2 = val2 * val2;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = tape_size(res2, val1, val2, val3);
    static_assert(size2 == 4);
}

TEST(adhoc2, tapeskip) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = valsum * val3;
    auto res = valsum * valprod;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 3);
}

TEST(adhoc2, tapecomplexd2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2 * adouble(2) + val3;
    auto valprod = val1 * val2 + adouble(3);
    auto valprod2 = val2 * val2 * adouble(4);
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = tape_size(res2, val1, val2, val3);
    static_assert(size2 == 4);
}

TEST(adhoc2, tapetrip) {
    adouble val1(1.);
    auto m1 = val1 * adouble(2.0);
    auto m2 = val1 * adouble(3.0);
    auto valprod = m1 * m2;

    constexpr std::size_t size = tape_size(valprod, val1);
    static_assert(size == 2);
}

TEST(adhoc2, tapetrip2) {
    adouble val1(1.);
    adouble val2(2.);
    auto m = val1 * val2;
    auto valprod = m * m;

    constexpr std::size_t size = tape_size(valprod, val1, val2);
    static_assert(size == 2);
}

TEST(adhoc2, tapecomplexdexp) {
    adouble val1(1.);
    adouble val2(2.);
    auto valsum = exp(val2);

    std::size_t size = tape_size(valsum, val2);
}

TEST(adhoc2, tapecomplexdexp2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + exp(val2) + val3;
    auto valprod = val1 * val2 + adouble(3);
    auto valprod2 = val2 * val2 * adouble(4);
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = tape_size(res2, val1, val2, val3);
    static_assert(size2 == 4);
}

// inline auto tuple_new() {
//     return std::make_tuple(adhoc<ID>(1.0), adhoc<ID>(2.0));
// }

// inline auto tuple_new2() {
//     return std::make_tuple(adhoc<ID>(1.0), adhoc<ID>(2.0));
// }

// void f(double) { std::cout << "1 arg" << std::endl; }
// void f(double, double) { std::cout << "2 arg" << std::endl; }
// void f(double, double, double) { std::cout << "3 arg" << std::endl; }
// void f(double, double, double, double) { std::cout << "4 arg" << std::endl; }
// void f(double, double, double, double, double) {
//     std::cout << "5 arg" << std::endl;
// }
// void f(double, double, double, double, double, double) {
//     std::cout << "6 arg" << std::endl;
// }

// #define DEBUG(...) __VA_OPT__(, ) __VA_ARGS__

// #define OP (
// #define CP )

// #define F(...) f OP 0 __VA_OPT__(, ) __VA_ARGS__ CP
// // #define F2(...) \
// //     std::make_tuple(adhoc<ID>(1.0 __VA_OPT__(CP, adhoc<ID> OP) __VA_ARGS__
// )) #define F2(A, B) std::make_tuple(adhoc<ID>(A CP, adhoc<ID> OP B))

// // ), adhoc<ID>(
// // CP, adhoc<ID>OP
// // std::make_tuple(adhoc<ID>(
// // std::make_tuple OP adhoc<ID> OP
// TEST(adhoc2, manyID) {
//     // constexpr std::size_t N = 10;
//     // std::array<double, N> in{};
//     // in[0] = 1.;
//     // in[1] = 2.;
//     // in[2] = 3.;
//     // in[3] = 4.;
//     // in[4] = 5.;
//     // in[5] = 6.;
//     // in[6] = 7.;
//     // in[7] = 8.;
//     // in[8] = 9.;
//     // in[9] = 10.;
//     // std::string test = "DEBUG(foo, bar, baz)";
//     // std::cout << test << std::endl;
//     // auto vars = tuple_new();
//     // std::cout << std::get<0>(vars).v() << std::endl;
//     // std::cout << std::get<1>(vars).v() << std::endl;
//     std::vector<double> v(10);
//     F();
//     F(v[0]);
//     F(v[0], v[1]);
//     F(v[0], v[1], v[2]);
//     F(v[0], v[1], v[3], v[4]);
//     F(v[0], v[1], v[3], v[4], v[5]);

//     auto tempvar = F2(1.0, 2.0);
//     std::cout << std::get<0>(tempvar).v();
//     std::cout << std::get<1>(tempvar).v();
//     auto tempvar2 = std::make_tuple(adhoc<ID>(), adhoc<ID>());
//     // F(v[0], v[1], v[3], v[4], v[5], v[6]);
//     // F(v[0], v[1], v[3], v[4], v[5], v[6], v[7]);
// }

} // namespace adhoc2