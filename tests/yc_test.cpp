#include <gtest/gtest.h>

#include <base.hpp>
#include <init.hpp>
#include <tape.hpp>

#include "call_price.hpp"
#include "type_name.hpp"
#include <constants_type.hpp>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
// #include <execution>
#include <iostream>
#include <tuple>
#include <type_traits>

#include <numeric>

namespace adhoc {

struct mycurve {
    std::vector<std::pair<int, double>> m_values;

    auto getDiscountFactorOld(const int &in) -> double {
        auto it = std::min(
            std::upper_bound(
                m_values.begin(), m_values.end(),
                std::make_pair(in, std::numeric_limits<double>::max())),
            m_values.end() - 1);

        double diff = static_cast<double>(in - (it - 1)->first) / 365.25;

        double df1 = (it - 1)->second;
        double df2 = it->second;
        double df = df2 / df1;
        double diff_total =
            static_cast<double>(it->first - (it - 1)->first) / 365.25;
        double rate = -std::log(df) / diff_total;
        return (it - 1)->second * std::exp(-rate * diff);
    }

    auto getDiscountFactor(const int &in) -> double {
        auto it = std::min(
            std::upper_bound(
                m_values.begin(), m_values.end(),
                std::make_pair(in, std::numeric_limits<double>::max())),
            m_values.end() - 1);

        // auto it = std::max(
        //     std::min(
        //         std::upper_bound(
        //             m_values.begin(), m_values.end(),
        //             std::make_pair(in, std::numeric_limits<double>::max())),
        //         m_values.end() - 1),
        //     m_values.begin() + 1);

        double df1 = (it - 1)->second;
        double df2 = it->second;
        double df = df2 / df1;
        double proportion = static_cast<double>(in - (it - 1)->first) /
                            static_cast<double>(it->first - (it - 1)->first);

        return df1 * std::exp(std::log(df) * proportion);
    }
};

enum class Extrapolation { Flat, Linear, LinearExponential };
enum class Interpolation { Flat, Linear, LinearExponential };

template <Extrapolation Left, Interpolation Middle, Extrapolation Right>
struct PiecewiseFuction {

    explicit PiecewiseFuction(std::vector<std::pair<double, double>> &values)
        : m_values(values) {
        m_rates.resize(m_values.size());
        m_rates.front() = 0;
        for (std::size_t i = 1; i < m_values.size(); i++) {
            m_rates[i] = std::log(m_values[i].second / m_values[i - 1].second) /
                         (m_values[i].first - m_values[i - 1].first);
        }
    }
    std::vector<std::pair<double, double>> m_values;

    std::vector<double> m_rates;

    auto getDiscountFactor(const double &in) -> double {
        auto it_middle = std::upper_bound(
            m_values.begin(), m_values.end(),
            std::make_pair(in, std::numeric_limits<double>::max()));

        auto it_right = std::min(it_middle, m_values.end() - 1);
        auto it_left = std::max(it_right - 1, m_values.begin());

        double df1 = it_left->second;

        long int d = std::distance(m_values.begin(), it_right);
        double proportion = (in - it_left->first);

        return df1 *
               std::exp(m_rates[static_cast<unsigned int>(d)] * proportion);
    }

    auto getDiscountFactor2(const double &in) -> double {
        auto it_middle = std::upper_bound(
            m_values.begin(), m_values.end(),
            std::make_pair(in, std::numeric_limits<double>::max()));

        auto it_right = std::min(it_middle, m_values.end() - 1);
        auto it_left = std::max(it_right - 1, m_values.begin());

        double df1 = it_left->second;

        long int d = std::distance(m_values.begin(), it_right);
        // double proportion = static_cast<double>(in - it_left->first);

        return df1 * std::exp(m_rates[static_cast<unsigned int>(d)] *
                              (static_cast<double>(in) -
                               static_cast<double>(it_left->first)));
    }
};

TEST(YC, CurveReplication) {
    // int start(2019, 1, 1);
    // bsl::vector<bdlt::Date> dates{ bdlt::Date(2020, 1, 1), bdlt::Date(2021,
    // 1, 1) }; icsu::DblVector dfs{ 0.9, 0.8 }; auto curve =
    // icpsc::StdPointSpotCurve::makeStd(start,
    //                                                icsk::EmbeddingTime::ET_ACT_365_25,
    //                                                ifte::InterpMethod::Type::STEP_FORWARD,
    //                                                dates,
    //                                                dfs,
    //                                                icsu::DblVector{},
    //                                                "");

    mycurve temp;
    temp.m_values.emplace_back(0, 1.0);
    // temp.m_values.emplace_back(1, 1.0);
    temp.m_values.emplace_back(1 + 365, 1.0);
    temp.m_values.emplace_back(1 + 2 * 365, 0.9);
    temp.m_values.emplace_back(1 + 3 * 365, 0.8);

    std::vector<std::pair<double, double>> values;
    // temp2.m_values.emplace_back(0, 1.0);
    // values.emplace_back(1 + 364, 1.0);
    values.emplace_back(1 + 365, 1.0);
    values.emplace_back(1 + 2 * 365, 0.9);
    values.emplace_back(1 + 3 * 365, 0.8);

    PiecewiseFuction<Extrapolation::Flat, Interpolation::LinearExponential,
                     Extrapolation::Linear>
        temp2(values);

    for (int i = 2; i < 1000; i++) {
        double df1 = temp2.getDiscountFactor(i);
        // double df1 = curve->getDiscountFactor(start + i);
        double df2 = temp.getDiscountFactor(i);
        EXPECT_NEAR(df1, df2, 1e-14);
        // std::cout << df1 << std::endl;
    }
}

struct LinearInterpolation {
    explicit LinearInterpolation(std::vector<std::pair<double, double>> &values)
        : m_values(values) {}
    std::vector<std::pair<double, double>> m_values;

    auto getValue(const double in) const -> double {
        auto it_middle =
            std::ranges::upper_bound(std::cbegin(m_values), std::cend(m_values),
                                     in, {}, &std::pair<double, double>::first);

        auto it_right =
            std::max(std::min(it_middle, std::prev(std::cend(m_values))),
                     std::next(std::cbegin(m_values)));
        auto it_left = it_right - 1;

        double slope = (it_right->second - it_left->second) /
                       (it_right->first - it_left->first);

        double increase = (in - it_left->first) * slope;
        return it_left->second + increase;
    }

    auto getValue2(const double in) const -> double {
        auto it_middle =
            std::ranges::upper_bound(std::cbegin(m_values), std::cend(m_values),
                                     in, {}, &std::pair<double, double>::first);

        auto it_right =
            std::max(std::min(it_middle, std::prev(std::cend(m_values))),
                     std::next(std::cbegin(m_values)));
        auto it_left = std::prev(it_right);

        auto pos = (in - it_left->first) / (it_right->first - it_left->first);
        return std::lerp(it_left->second, it_right->second, pos);
    }
};

auto getValueStandalone(const double in,
                        std::vector<std::pair<double, double>> const &m_values)
    -> double {
    auto it_middle =
        std::ranges::upper_bound(std::cbegin(m_values), std::cend(m_values), in,
                                 {}, &std::pair<double, double>::first);

    auto it_right =
        std::max(std::min(it_middle, std::prev(std::cend(m_values))),
                 std::next(std::cbegin(m_values)));
    auto it_left = std::prev(it_right);

    auto pos = (in - it_left->first) / (it_right->first - it_left->first);
    return std::lerp(it_left->second, it_right->second, pos);
}

struct LinearInterpolation2 {
    explicit LinearInterpolation2(
        std::vector<std::pair<double, double>> &values) {
        m_values.resize(values.size());
        m_values.front() =
            std::make_tuple(values[0].first, values[0].second, 0);

        for (std::size_t i = 1; i < m_values.size(); i++) {
            m_values[i] =
                std::make_tuple(values[i].first, values[i].second,
                                1.0 / (values[i].first - values[i - 1].first));
        }
    }
    std::vector<std::tuple<double, double, double>> m_values;

    auto getValue(const double in) -> double {
        auto it_middle = std::ranges::upper_bound(
            std::cbegin(m_values), std::cend(m_values), in, {},
            [](const auto &t) { return std::get<0>(t); });

        auto it_right =
            std::max(std::min(it_middle, std::prev(std::cend(m_values))),
                     std::next(std::cbegin(m_values)));
        auto i = std::distance(std::cbegin(m_values), it_right);

        auto pos =
            (in - std::get<0>(m_values[static_cast<std::size_t>(i - 1)])) *
            std::get<2>(m_values[static_cast<std::size_t>(i)]);
        return std::lerp(std::get<1>(m_values[static_cast<std::size_t>(i - 1)]),
                         std::get<1>(m_values[static_cast<std::size_t>(i)]),
                         pos);
    }

    auto getValue2(const double in) const -> double {
        auto it_middle = std::ranges::upper_bound(
            std::cbegin(m_values), std::cend(m_values), in, {},
            [](const auto &t) { return std::get<0>(t); });

        auto it_right =
            std::max(std::min(it_middle, std::prev(std::cend(m_values))),
                     std::next(std::cbegin(m_values)));
        auto it_left = std::prev(it_right);

        auto pos = (in - std::get<0>(*it_left)) * std::get<2>(*it_right);
        return std::lerp(std::get<1>(*it_left), std::get<1>(*it_right), pos);
    }
};

struct LinearInterpolation3 {
    explicit LinearInterpolation3(
        std::vector<std::pair<double, double>> &values) {
        m_values.resize(values.size());
        m_values.front() =
            std::make_tuple(values[0].first, values[0].second, 0);

        for (std::size_t i = 1; i < m_values.size(); i++) {
            m_values[i] =
                std::make_tuple(values[i].first, values[i].second,
                                (values[i].second - values[i - 1].second) /
                                    (values[i].first - values[i - 1].first));
        }
    }
    std::vector<std::tuple<double, double, double>> m_values;

    auto getValue(const double in) -> double {
        auto it_middle = std::ranges::upper_bound(
            std::cbegin(m_values), std::cend(m_values), in, {},
            [](const auto &t) { return std::get<0>(t); });

        auto it_right =
            std::max(std::min(it_middle, std::prev(std::cend(m_values))),
                     std::next(std::cbegin(m_values)));
        auto it_left = std::prev(it_right);

        // double slope = (it_right->second - it_left->second) /
        //                (it_right->first - it_left->first);

        double increase = (in - std::get<0>(*it_left)) * std::get<2>(*it_right);
        return std::get<1>(*it_left) + increase;
    }

    // auto getValue2(const double in) const -> double {
    //     auto it_middle = std::upper_bound(
    //         std::cbegin(m_values), std::cend(m_values),
    //         std::make_tuple(in, std::numeric_limits<double>::max(),
    //                         std::numeric_limits<double>::max()));

    //     auto it_right =
    //         std::max(std::min(it_middle, std::prev(std::cend(m_values))),
    //                  std::next(std::cbegin(m_values)));
    //     auto it_left = std::prev(it_right);

    //     auto pos = (in - std::get<0>(*it_left)) * std::get<2>(*it_right);
    //     return std::lerp(std::get<1>(*it_left), std::get<1>(*it_right), pos);
    // }
};

TEST(YC, LinearInterpolation) {
    std::vector<std::pair<double, double>> values;
    values.emplace_back(0, 1.0);
    values.emplace_back(1, 3.0);
    values.emplace_back(5, 2);
    values.emplace_back(10, 0);
    LinearInterpolation temp(values);

    double start = -1.0;
    double increment = 0.1;

    std::vector<double> results = {
        -1, -0.8,  -0.6, -0.4,  -0.2, 0,     0.2,  0.4,   0.6, 0.8,
        1,  1.2,   1.4,  1.6,   1.8,  2,     2.2,  2.4,   2.6, 2.8,
        3,  2.975, 2.95, 2.925, 2.9,  2.875, 2.85, 2.825, 2.8, 2.775};
    // double val = start;

    LinearInterpolation2 temp2(values);
    for (int i = 0; i < 30; i++) {
        double val = start + i * increment;
        double res = temp.getValue(val);
        double res2 = temp.getValue2(val);
        double res3 = temp2.getValue(val);
        EXPECT_NEAR(res, results[static_cast<std::size_t>(i)], 1e-14);
        EXPECT_NEAR(res2, results[static_cast<std::size_t>(i)], 1e-14);
        EXPECT_NEAR(res3, results[static_cast<std::size_t>(i)], 1e-14);
        val += increment;
    }
}

TEST(YC, LinearInterpolationTime) {
    std::vector<std::pair<double, double>> values;
    values.emplace_back(0, 1.0);
    values.emplace_back(1, 3.0);
    values.emplace_back(5, 2);
    values.emplace_back(10, 0);

    double start = -1.0;
    double increment = 0.0001;
    int iters = 100000;
    double sum = 0;
    auto starttime = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    // {
    //     starttime = std::chrono::steady_clock::now();
    //     sum = 0;
    //     LinearInterpolation temp(values);
    //     for (int i = 0; i < iters; i++) {
    //         double val = start + i * increment;
    //         double res = temp.getValue(val);
    //         sum += res;
    //     }
    //     end = std::chrono::steady_clock::now();

    //     std::cout << sum << std::endl;
    //     std::cout << "Elapsed time in microseconds: "
    //               << std::chrono::duration_cast<std::chrono::microseconds>(
    //                      end - starttime)
    //                      .count()
    //               << " µs" << std::endl;
    // }

    {
        starttime = std::chrono::steady_clock::now();
        sum = 0;
        LinearInterpolation temp2(values);
        for (int i = 0; i < iters; i++) {
            double val = start + i * increment;
            double res = temp2.getValue2(val);
            sum += res;
        }
        end = std::chrono::steady_clock::now();

        std::cout << sum << std::endl;
        std::cout << "Elapsed time in microseconds: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
    }

    {
        starttime = std::chrono::steady_clock::now();
        sum = 0;
        LinearInterpolation2 temp3(values);
        for (int i = 0; i < iters; i++) {
            double val = start + i * increment;
            double res = temp3.getValue(val);
            sum += res;
        }
        end = std::chrono::steady_clock::now();

        std::cout << sum << std::endl;
        std::cout << "Elapsed time in microseconds: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
    }

    {
        starttime = std::chrono::steady_clock::now();
        sum = 0;
        LinearInterpolation2 temp4(values);
        for (int i = 0; i < iters; i++) {
            double val = start + i * increment;
            double res = temp4.getValue2(val);
            sum += res;
        }
        end = std::chrono::steady_clock::now();

        std::cout << sum << std::endl;
        std::cout << "Elapsed time in microseconds: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
    }

    {
        starttime = std::chrono::steady_clock::now();
        sum = 0;
        LinearInterpolation3 temp4(values);
        for (int i = 0; i < iters; i++) {
            double val = start + i * increment;
            double res = temp4.getValue(val);
            sum += res;
        }
        end = std::chrono::steady_clock::now();

        std::cout << sum << std::endl;
        std::cout << "Elapsed time in microseconds: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
    }

    {
        starttime = std::chrono::steady_clock::now();
        sum = 0;
        LinearInterpolation temp(values);
        for (int i = 0; i < iters; i++) {
            double val = start + i * increment;
            double res = temp.getValue(val);
            sum += res;
        }
        end = std::chrono::steady_clock::now();

        std::cout << sum << std::endl;
        std::cout << "Elapsed time in microseconds: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
    }

    {
        starttime = std::chrono::steady_clock::now();
        sum = 0;
        // LinearInterpolation temp(values);
        for (int i = 0; i < iters; i++) {
            double val = start + i * increment;
            double res = getValueStandalone(val, values);
            sum += res;
        }
        end = std::chrono::steady_clock::now();

        std::cout << sum << std::endl;
        std::cout << "Elapsed time in microseconds: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
    }
}

template <typename T> class vector;

template <class Loc, class I> struct forward_iterator {
    // auto operator->() -> I* { return &I{}; }
    // auto operator->() const -> I* { return I{}; }
    auto operator*() -> I { return I{}; }
    auto operator*() const -> I { return I{}; }
};

class begin_it_t;
class end_it_t;
// template <class T> auto begin(vector<T> /* a */) {
//     return forward_iterator<T>();
// }

// template <class T> auto end(vector<T> /* a */) { return
// forward_iterator<T>(); }

// template <class T> auto cbegin(vector<T> /* a */) {
//     return forward_iterator<T>();
// }

// template <class T> auto cend(vector<T> /* a */) {
//     return forward_iterator<T>();
// }

// template <std::size_t N> struct double_t : public Base<double_t<N>> {};

template <class V> class begin_t : public forward_iterator<begin_it_t, V> {};
template <class V> class end_t : public forward_iterator<end_it_t, V> {};

template <typename T> class vector {
  public:
    explicit vector() = default;
    [[nodiscard]] auto begin() const -> begin_t<T> { return begin_t<T>(); }
    [[nodiscard]] auto end() const -> end_t<T> { return end_t<T>(); }
};

template <class V> auto begin(V /* in */) { return begin_t<V>(); }
template <class V> auto end(V /* in */) { return end_t<V>(); }

// template <std::size_t N> vector<N>::vector() {}

template <class V, class Init, class Lambda>
class accumulate_t : public Base<accumulate_t<V, Init, Lambda>> {
  public:
    // explicit accumulate_t();
    void print();
    [[nodiscard]] auto v() const noexcept -> double { return 0; }
};

template <class V, class Init, class Lambda>
void accumulate_t<V, Init, Lambda>::print() {
    std::cout << "I am in" << std::endl;
}

// we only support begin() and end() as idiom
template <typename VT, class T, class BinaryOperation>
auto accumulate(begin_t<vector<VT>> /* first */, end_t<vector<VT>> /* last */,
                T init, BinaryOperation op) {
    return accumulate_t<vector<VT>, T, decltype(op(init, VT{}))>();
}

template <typename VT, class T, class BinaryOperation>
auto reduce(begin_t<vector<VT>> /* first */, end_t<vector<VT>> /* last */,
            T init, BinaryOperation op) {
    return accumulate_t<vector<VT>, T, decltype(op(init, VT{}))>();
}

// template <class T, std::size_t N>
// struct vector : public adhoc::Base<vector<T, N>> {};

// template <class V, class Init, class Lambda>
// void accumulate_t<V, Init, Lambda>::print() {
//     std::cout << "I am in" << std::endl;
// }

// // we only support begin() and end() as idiom
// template <std::size_t N, class T, class BinaryOperation>
// auto accumulate(begin_t<vector<N>> /* first */, end_t<vector<N>> /* last */,
//                 T init, BinaryOperation op) {
//     return accumulate_t<vector<N>, T, decltype(op(init, double_t<N>{}))>();
// }
namespace std2 {
template <class Container, class UnaryOperation>
auto transform_container(const Container &ic, UnaryOperation f) -> Container {
    Container oc;
    std::transform(std::begin(ic), std::end(ic), std::inserter(oc, oc.end()),
                   f);
    return oc;
}
} // namespace std2

template <class Container, class UnaryOperation> class transform_container_t {};
template <class Container, class UnaryOperation>
auto transform_container(const Container & /* ic */, UnaryOperation /* f */) {
    return transform_container_t<Container, UnaryOperation>{};
}

TEST(YC, MapReduce) {
    vector<double_t<5>> data_adhoc;
    // vector<double> data2_adhoc;

    auto lambdafunc2 = []<typename T>(T const &lhs) { return lhs * lhs; };
    auto data2_adhoc = transform_container(data_adhoc, lambdafunc2);

    std::vector<double> data(10);
    std::iota(std::begin(data), std::end(data), 0.);

    // using std2::transform_container;
    auto reducfunc = [](double const &lhs, double const &rhs) {
        return lhs + rhs;
    };
    auto lambdafunc = [](double const &lhs) { return lhs * lhs; };
    double res = std::transform_reduce(std::cbegin(data), std::cend(data), 0,
                                       reducfunc, lambdafunc);
    std::cout << res << std::endl;
}

template <class V, class Init, class BinaryReductionOp, class UnaryTransformOp>
class transform_reduce_t
    : public Base<
          transform_reduce_t<V, Init, BinaryReductionOp, UnaryTransformOp>> {
  public:
    void print();
    [[nodiscard]] auto v() const noexcept -> double { return 0; }
};

template <class V, class Init, class BinaryReductionOp, class UnaryTransformOp>
void transform_reduce_t<V, Init, BinaryReductionOp, UnaryTransformOp>::print() {
    std::cout << "I am in" << std::endl;
}

template <typename VT, class T, class BinaryReductionOp, class UnaryTransformOp>
auto transform_reduce2(begin_t<VT> /* first */, end_t<VT> /* last */, T init,
                       BinaryReductionOp reduce, UnaryTransformOp transform) {
    return transform_reduce_t<vector<VT>, T,
                              decltype(reduce(init, transform(VT{}))),
                              decltype(transform(VT{}))>();
}

template <class A, class B, class T>
class lerp_t : public Base<lerp_t<A, B, T>> {};

template <class A, class B, class T>
auto lerp(A /* a */, B /* b */, T /* t */) {
    return lerp_t<A, B, T>();
}

template <typename T1, typename T2> struct pair {
    pair() = default;
    T1 first;
    T2 second;
};

// template <class I, class T, class Proj>
// class upper_bound_t : public Base<upper_bound_t<I, T, Proj>> {};

template <class I, class T, class Proj>
class upper_bound_t : public forward_iterator<upper_bound_t<I, T, Proj>, I> {};

template <class I, class T, class Comp = std::pair<double, double>,
          class Proj = std::pair<double, double>>
auto upper_bound(begin_t<I> /* first */, end_t<I> /* last */,
                 const T & /* value */, Comp /* comp */ = {}, Proj proj = {}) {
    // I val{};
    return forward_iterator<upper_bound_t<I, T, decltype(proj(I{}))>, I>{};
};

TEST(YC, MapReduceBond) {
    std::vector<std::pair<double, double>> payments = {
        {1., 0.02}, {2., 0.02}, {3., 0.01}, {4., 0.01}, {5., 0.03}, {6., 1.03}};

    std::vector<std::pair<double, double>> values;
    values.emplace_back(0, 1.0);
    values.emplace_back(1, 0.98);
    values.emplace_back(5, 0.97);
    values.emplace_back(10, 0.94);

    using std::begin;
    using std::cbegin;
    using std::cend;
    using std::end;
    using std::transform_reduce;
    auto LinearInterp = [&values]<typename T>(const T in) {
        using std::ranges::upper_bound;
        using std::next;
        using std::prev;
        using std::lerp;
        using std::min;
        using std::max;
        using std::get;
        auto it_middle = upper_bound(begin(values), end(values), in, {},
                                     [](const auto &t) { return get<0>(t); });

        auto it_right =
            max(min(it_middle, prev(end(values))), next(begin(values)));
        auto it_left = prev(it_right);

        auto pos =
            (in - get<0>(*it_left)) / (get<0>(*it_right) - get<0>(*it_left));
        return lerp(get<1>(*it_left), get<1>(*it_right), pos);
    };

    auto reducfunc = []<typename T1, typename T2>(
                         T1 const &lhs, T2 const &rhs) { return lhs + rhs; };

    auto lambdafunc = [&]<typename T>(T const &lhs) {
        return lhs.second * LinearInterp(std::get<0>(lhs));
    };

    //
    // using constants::CD;
    // using constants::encode;
    // CD<encode(0)> init;
    // double initd = init;
    double res = transform_reduce(cbegin(payments), cend(payments), 0.,
                                  reducfunc, lambdafunc);

    EXPECT_NEAR(res, 1.080645, 1e-14);

    using constants::CD;
    using constants::encode;
    CD<encode(0)> init;
    CD<encode(3)> value;
    vector<pair<double_t<3>, double_t<4>>> values2;
    auto itbegin = values2.begin();
    auto itend = values2.end();
    double_t<10> in;

    // auto reducfunc = []<typename T1, typename T2>(
    //                      T1 const &lhs, T2 const &rhs) { return lhs + rhs; };

    auto projfunc = []<typename T>(T const &t) { return t.first; };

    auto its = upper_bound(itbegin, itend, in, {}, projfunc);
    auto val = (*its).first;
    std::cout << type_name<decltype(val)>() << std::endl;
    // auto val = its->first;
    // its->

    adhoc::pair<double_t<3>, double_t<4>> singlevalue2{};
    auto LinearInterpSimple = [&value, &singlevalue2]<typename T>(const T in) {
        return lerp(singlevalue2.first, singlevalue2.second, in);
    };

    auto simpletransform = [&LinearInterpSimple]<typename T>(T const &lhs) {
        return lhs.first * LinearInterpSimple(lhs.second);
    };

    auto simplereduc = []<typename T1, typename T2>(
                           T1 const &lhs, T2 const &rhs) { return lhs + rhs; };
    vector<pair<double_t<5>, double_t<6>>> payments2;
    auto res2 = adhoc::transform_reduce2(cbegin(payments2), cend(payments2),
                                         init, simplereduc, simpletransform);

    // std::cout << type_name<decltype(res2)>() << std::endl;
}

} // namespace adhoc
