#include <adhoc2.hpp>
#include <evaluate.hpp>
#include <tape_size.hpp>
#include <utils.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, Initial) {
    adouble val(1.);
    auto valplus = val + adouble(1.0);
    auto const valmult = valplus * adouble(2.0);
    auto valmult2 = (valmult + adouble(2.0)) * adouble(3.0);
}

TEST(adhoc2, depends2) {
    adouble val1(1.);
    adouble val2(2.);
    constexpr auto res =
        decltype(val1)::depends3<decltype(val1), decltype(val1),
                                 decltype(val2)>();
    EXPECT_EQ(res, 2);
}

TEST(adhoc2, complexdepends2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = val1 * val2;
    auto t1 = valsum * valprod;
    auto t2 = valprod + t1;

    constexpr auto res = decltype(t2)::depends3<decltype(val1), decltype(val2),
                                                decltype(val3)>();
    EXPECT_EQ(res, 2);
}

TEST(adhoc2, complexdepends3) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = val1 * val2;
    auto t1 = valsum * valprod;
    auto t2 = valprod + t1;

    constexpr auto res = decltype(t2)::depends<decltype(t1)>();
    EXPECT_EQ(res, true);
}

template <int N> constexpr auto first_n_fibs() -> std::array<double, N> {
    std::array<double, N> ret{};
    for (std::size_t i = 0; i < N; i++)
        ret[i] = 1.0;
    return ret;
}

TEST(adhoc2, complexd2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = val1 * val2;
    auto t1 = valsum * valprod;
    auto t2 = valprod + t1;

    constexpr auto a = first_n_fibs<3>();
    // constexpr auto res = t2.d2(val1, val2);
    EXPECT_EQ(a.size(), 3);
}

template <std::size_t N> class sometype {};

TEST(adhoc2, position) {
    double temp = 0;
    constexpr auto res1 = has_type<int, float, double, int, unsigned int>(temp);
    constexpr auto res2 = has_type<int, float, int, unsigned int>(temp);
    static_assert(res1);
    static_assert(!res2);
    constexpr auto res3 = idx_type<int, float, double, int, unsigned int>(temp);
    static_assert(res3 == 2);

    constexpr auto res4 = has_type2<double>();
    static_assert(!res4);
}

TEST(adhoc2, derivcomplexdexp2) {
    adouble x(1.);
    adouble y(2.);
    adouble z(2.);
    auto expy = exp(y);
    auto s1 = x + expy;
    auto s2 = s1 + z;
    auto xy = x * y;
    auto xyp3 = xy + adouble(3);
    std::cout << xyp3.input2().v() << std::endl;
    auto yy = y * y;
    auto yyt4 = yy + adouble(4);
    auto p = xyp3 * yyt4;
    auto res = s2 * p;

    constexpr std::size_t size = tape_size(res, x, y, z);
    static_assert(size == 5);
    std::array<double, size> vals = {};
    std::cout << s2.v() << std::endl;
    std::cout << yyt4.v() << std::endl;

    vals[0] = res.d1(); // valsum
    vals[1] = res.d2(); // valprod3

    vals[2] = vals[0];
    vals[0] *= s2.d1(); // tmp1
    vals[2] *= s2.d2(); // val3

    vals[3] = vals[0];
    vals[0] *= s1.d1(); // val1
    vals[3] *= s1.d2(); // exp2

    vals[3] *= expy.d(); // val2

    vals[4] = vals[1];
    vals[1] *= p.d1(); // valprod
    vals[4] *= p.d2(); // valprod2

    vals[1] *= xyp3.d1(); // val12

    vals[0] += vals[1] * xy.d1(); // val1
    vals[3] += vals[1] * xy.d2(); // val2

    vals[4] *= yyt4.d1(); // val22

    vals[3] += vals[1] * yy.d1(); // val2
    vals[3] += vals[1] * yy.d2(); // val2

    // EXPECT_NEAR(vals[0], 412.4497951657808, 1e-9);
    // EXPECT_NEAR(vals[1], 1588.4738734117946, 1e-9);
    // EXPECT_NEAR(vals[2], 80., 1e-9);

    std::cout << vals[0] << ", " << vals[3] << ", " << vals[2] << std::endl;
    std::cout << vals[0] << ", " << vals[1] << ", " << vals[2] << ", "
              << vals[3] << ", " << vals[4] << std::endl;
    auto res2 = p * s2;
    constexpr std::size_t size2 = tape_size(res2, x, y, z);
    static_assert(size2 == 4);
}

TEST(adhoc2, derivcomplexdexpevaluate) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto exp2 = exp(val2);
    auto tmp1 = val1 + exp2;
    auto valsum = tmp1 + val3;
    auto val12 = val1 * val2;
    auto valprod = val12 + adouble(3);
    auto val22 = val2 * val2;
    auto valprod2 = val22 * adouble(4);
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    auto derivatives = evaluate(res, val1, val2, val3);
    static_assert(derivatives.size() == 3);
    std::cout << derivatives[0] << std::endl;
    std::cout << derivatives[1] << std::endl;
    std::cout << derivatives[2] << std::endl;
}

// TEST(adhoc2, derivcomplexdexpevaluate) {
//     adouble val1(1.);
//     adouble val2(2.);
//     auto valprod3 = (val1 * val2) + (val1 * val1);

//     auto derivatives = evaluate(valprod3, val1, val2);
// }

TEST(adhoc2, derivcomplexdexpevaluatesimple) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod3 = (val1 * val2) + (val1);

    auto derivatives = evaluate(valprod3, val1, val2);
}

TEST(adhoc2, derivativezero) {
    adouble val1(1.);
    adouble val2(1.);
    auto derivatives = evaluate(val1, val2);
    static_assert(derivatives.size() == 0);
}

TEST(adhoc2, derivativeexp) {
    adouble val1(1.);
    auto temp = exp(val1);
    auto derivatives = evaluate(temp, val1);
    static_assert(derivatives.size() == 1);
}

TEST(adhoc2, derivativecos) {
    adouble val1(1.);
    auto temp = cos(val1);
    auto derivatives = evaluate(temp, val1);
    static_assert(derivatives.size() == 1);
}

TEST(adhoc2, derivativeadd) {
    adouble val1(1.);
    adouble val2(1.);
    auto temp = val1 + val2;
    auto derivatives = evaluate(temp, val1, val2);
    static_assert(derivatives.size() == 2);
}

TEST(adhoc2, derivativemul) {
    adouble val1(1.);
    adouble val2(1.);
    auto temp = val1 * val2;
    auto derivatives = evaluate(temp, val1, val2);
    static_assert(derivatives.size() == 2);
}

TEST(adhoc2, derivative) {
    adouble val1(1.);
    auto valexp = exp(val1);
    auto valcos = cos(valexp);
    constexpr std::size_t size = tape_size(valcos, val1);
    static_assert(size == 1);

    double derivative = valcos.d();
    derivative *= valexp.d();
    // www.wolframalpha.com/input?i=d%2Fdx+cos%28exp%28x%29%29+%7C+x%3D1
    EXPECT_NEAR(derivative, -1.116619317445013, 1e-10);
}

TEST(adhoc2, DerivativeUniv) {
    adouble x(1.);
    auto result = evaluate(cos(exp(x)), x);

    static_assert(result.size() == 1);
    // from sympy import *
    // x = Symbol('x')
    // f = cos(exp(x))
    // dfdx = f.diff(x)
    // dfdx = lambdify(x, dfdx)
    // print(dfdx(1))
    EXPECT_NEAR(result[0], -1.1166193174450132, 1e-10);
}

auto makeTemp() {
    adouble x(1.);
    auto valexp = exp(x);
    auto valcos = cos(valexp);
    return std::tuple<decltype(valcos), decltype(x)>{valcos, x};
}

TEST(adhoc2, DerivativeUnivScope) {
    // this shows why we need deep copies. thisexample would not work on release
    // if we would use references
    auto temp = makeTemp();
    auto result = evaluate(std::get<0>(temp), std::get<1>(temp));
    static_assert(result.size() == 1);
    EXPECT_NEAR(result[0], -1.1166193174450132, 1e-10);
}

TEST(adhoc2, DerivativeMult) {
    adouble x(1.5);
    adouble y(2.5);
    auto result = evaluate(x * y, x, y);

    static_assert(result.size() == 2);
    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = x * y
    // dfdx = f.diff(x)
    // dfdy = f.diff(y)
    // dfdx = lambdify([x, y], dfdx)
    // dfdy = lambdify([x, y], dfdy)
    // print(dfdx(1.5, 2.5))
    // print(dfdy(1.5, 2.5))
    EXPECT_NEAR(result[0], 2.5, 1e-10);
    EXPECT_NEAR(result[1], 1.5, 1e-10);
    std::cout << result[0] << std::endl;
    std::cout << result[1] << std::endl;
}

TEST(adhoc2, DerivativeMultSame) {
    adouble x(1.5);
    auto result = evaluate(x * x, x);

    static_assert(result.size() == 1);
    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = x * y
    // dfdx = f.diff(x)
    // dfdy = f.diff(y)
    // dfdx = lambdify([x, y], dfdx)
    // dfdy = lambdify([x, y], dfdy)
    // print(dfdx(1.5, 2.5))
    // print(dfdy(1.5, 2.5))
    EXPECT_NEAR(result[0], 3.0, 1e-10);

    auto f = x * x;
    constexpr std::size_t size = tape_size(f, x);
    static_assert(size == 1);
}

// TEST(adhoc2, Derivative2nd) {
//     adouble x(2);
//     adouble y(3);
//     adouble z(5);
//     auto x2 = x * x;
//     auto y2 = y * y;
//     auto z2 = z * z;
//     auto temp = x2 * y2;
//     auto temp2 = temp * z2;
//     auto res = temp * temp2;

//     // f1(temp, temp2) = temp * temp2
//     // df1dtemp = dtemp/dtemp * temp2 + temp * dtemp2/dtemp
//     // df1dtemp_2 = dtemp/dtemp_2 * temp2 + temp * dtemp2/dtemp_2 +
//     dtemp/dtemp * dtemp2/dtemp double dresdtemp = temp2.v(); double
//     dresdtemp2 = temp.v(); double d2_res_d_temp_2 = 0; double
//     d2_res_d_temp2_2 = 0; double d2_res_dtemp2_dtemp = 1;

//     // f2(temp, z2) = temp * z2
//     // f2(temp, z2) = temp * (temp * z2) = f1(temp, f2(temp, z2))
//     // dfdtemp_2 = dtemp/dtemp_2 * (temp * z2) + temp * d(temp * z2)/dtemp_2
//     + dtemp/dtemp * d(temp * z2)/dtemp dresdtemp += dresdtemp2 * z2.v();
//     double d2resdz2 = dresdtemp2 * temp.v();
//     d2_res_d_temp_2 += 0;
//     double d2_res_d_z2_2 = 0;
//     double d2_res_dz2_dtemp = 0;
// }

} // namespace adhoc2