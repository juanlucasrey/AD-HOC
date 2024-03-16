#include <adhoc.hpp>
#include <constants_type.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

#include <random>

namespace adhoc {

auto transform_in_place_m(std::vector<double> &in, double val)
    -> std::vector<double> & {
    std::transform(in.begin(), in.end(), in.begin(), [val](auto &&PH1) {
        return std::multiplies<double>()(std::forward<decltype(PH1)>(PH1), val);
    });
    return in;
}

auto transform_in_place_a(std::vector<double> &in, double val)
    -> std::vector<double> & {
    std::transform(in.begin(), in.end(), in.begin(), [val](auto &&PH1) {
        return std::plus<double>()(std::forward<decltype(PH1)>(PH1), val);
    });

    return in;
}

// template <class I1, class I2, class I3>
// double lambdafunc(I1 &in, I2 mul = I2(), I3 add = I3()) {
//     return in * mul + add;
// }
template <class V> class begin_t {};
template <class V> class end_t {};

template <std::size_t N> class vector {
  public:
    explicit vector();
    [[nodiscard]] auto begin() const -> begin_t<vector<N>> {
        return begin_t<vector<N>>();
    }
    [[nodiscard]] auto end() const -> end_t<vector<N>> {
        return end_t<vector<N>>();
    }
};

template <std::size_t N> vector<N>::vector() {}

template <class V, class Init, class Lambda>
class accumulate_t : public Base<accumulate_t<V, Init, Lambda>> {
  public:
    // explicit accumulate_t();
    void print();
    [[nodiscard]] auto v() const noexcept -> double { return 0; }
};

// template <class Input1, class Input2>
// accumulate_t<Input1, Input2>::accumulate_t() {}

template <class V, class Init, class Lambda>
void accumulate_t<V, Init, Lambda>::print() {
    std::cout << "I am in" << std::endl;
}

// we only support begin() and end() as idiom
template <std::size_t N, class T, class BinaryOperation>
auto accumulate(begin_t<vector<N>> /* first */, end_t<vector<N>> /* last */,
                T init, BinaryOperation op) {
    return accumulate_t<vector<N>, T, decltype(op(init, double_t<N>{}))>();
}

#if __cplusplus >= 202002L
TEST(AD, BlackScholesSimulationadhoc) {
    double_t<0> S;
    double_t<1> K;
    double_t<2> v;
    double_t<3> T;
    vector<4> rndnmbrs;

    double_t<5> n_times;
    auto dt = T / n_times;

    auto adjvol = v * sqrt(dt);

    using constants::CD;
    using constants::encode;
    auto add = CD<encode(-0.5)>() * adjvol * adjvol;

    auto lambdafunc = [adjvol, add]<typename T1, typename T2>(T1 const &lhs,
                                                              T2 const &rhs) {
        return lhs + (rhs * adjvol + add);
    };

    CD<encode(0)> init;
    auto result =
        accumulate(rndnmbrs.begin(), rndnmbrs.end(), init, lambdafunc);
    std::cout << type_name<decltype(result)>() << std::endl;

    auto fwd = S * exp(result);
}

TEST(AD, BlackScholesSimulationdouble) {

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using namespace constants;

    double S = stock_distr(generator);
    double K = stock_distr(generator);
    double v = vol_distr(generator);
    double T = time_distr(generator);
    using std::accumulate;
    using std::exp;
    using std::sqrt;
    using std::vector;
    vector<double> rndnmbrs;

    double n_times(20);
    auto dt = T / n_times;

    auto adjvol = v * sqrt(dt);

    auto add = CD<encode(-0.5)>() * adjvol * adjvol;

    auto lambdafunc = [adjvol, add]<typename T1, typename T2>(T1 const &lhs,
                                                              T2 const &rhs) {
        return lhs + (rhs * adjvol + add);
    };

    double init(0);
    auto result =
        accumulate(rndnmbrs.begin(), rndnmbrs.end(), init, lambdafunc);

    auto fwd = S * exp(result);
    std::cout << fwd << std::endl;
}

TEST(AD, BlackScholesSimulation) {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    double S = stock_distr(generator);
    double K = stock_distr(generator);
    double v = vol_distr(generator);
    double T = time_distr(generator);

    std::size_t n_times = 20;
    double dt = T / static_cast<double>(n_times);

    std::mt19937 rng(123);
    std::normal_distribution<double> dd;

    std::size_t n_paths = 10000;

    // using namespace BloombergLP::qfd::quantmath::statistics;
    // accumulator<> acc(true);

    double adjvol = v * std::sqrt(dt);

    double mean = 0;
    for (std::size_t i = 0; i < n_paths; ++i) {
        // double fwd = S;
        // for (std::size_t j = 0; j < n_times; ++j) {
        //     fwd *= std::exp(dd(rng) * adjvol - 0.5 * adjvol * adjvol);
        // }

        std::vector<double> rndnmbrs(n_times);
        for (std::size_t j = 0; j < n_times; ++j) {
            rndnmbrs[j] = dd(rng);
        }

        double cumul_var = 0.;
        double mul = adjvol;
        double add = -0.5 * adjvol * adjvol;
        for (std::size_t j = 0; j < n_times; ++j) {
            cumul_var += rndnmbrs[j] * mul + add;
        }

        // auto lambda = [adjvol](auto && PH1) { return
        // std::multiplies<double>()(std::forward<decltype(PH1)>(PH1), adjvol);
        // }; std::cout << decltype(lambda) << std::endl;
        // auto newvec = transform_in_place_m(rndnmbrs, adjvol);
        // auto newvec2 = transform_in_place_a(newvec, -0.5 * adjvol * adjvol);
        // std::transform(
        //     rndnmbrs.begin(), rndnmbrs.end(), rndnmbrs.begin(),
        //     std::bind(std::multiplies<double>(), std::placeholders::_1, 3));

        auto lambdafunc = [mul, add]<typename T>(T const &lhs, T const &rhs) {
            return lhs + rhs * mul + add;
        };
        using std::accumulate;
        auto result =
            accumulate(rndnmbrs.begin(), rndnmbrs.end(), 0., lambdafunc);

        double fwd = S * std::exp(result);

        double price = fwd > K ? fwd - K : 0;

        mean += price;
    }

    mean /= static_cast<double>(n_paths);
    std::cout << mean << std::endl;

    // double adjvol = 1;
    auto lambda = [adjvol](auto &&PH1) {
        return std::multiplies<double>()(std::forward<decltype(PH1)>(PH1),
                                         adjvol);
    };
    std::cout << type_name<decltype(lambda)>() << std::endl;
    // std::cout << decltype(lambda) << std::endl;
}
#endif

TEST(AD, accumulate) {
    std::vector<double> test{0.1, 0.2, 0.3};
    auto result = std::accumulate(test.begin(), test.end(), 0);
    std::cout << result << std::endl;
}

} // namespace adhoc