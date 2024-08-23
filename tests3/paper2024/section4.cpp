#include <calc_tree.hpp>
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "../call_price.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <random>

namespace adhoc3 {

namespace {

template <class D> inline auto pdf_n(D x) {
    using std::exp;
    using namespace constants;
    constexpr double one_over_root_two = 0.39894228040143265;
    return CD<encode(one_over_root_two)>() * exp(CD<encode(-0.5)>() * x * x);
}

template <class I1, class I2, class I3, class I4>
inline auto vega(I1 S, I2 K, I3 v, I4 T) {
    using std::log;
    using std::sqrt;
    using namespace constants;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<encode(0.5)>();
    return S * pdf_n(d1) * sqrt(T);
}

inline auto vanna(double S, double K, double v, double T) {
    auto totalvol = v * std::sqrt(T);
    auto d1 = std::log(S / K) / totalvol + totalvol * 0.5;
    auto d2 = d1 - totalvol;
    return -pdf_n(d1) * d2 / v;
}

inline auto volga(double S, double K, double v, double T) {
    auto totalvol = v * std::sqrt(T);
    auto d1 = std::log(S / K) / totalvol + totalvol * 0.5;
    auto d2 = d1 - totalvol;
    return S * pdf_n(d1) * d1 * d2 * T / totalvol;
}

inline auto call_price_vega_vanna_volga(double S, double K, double v,
                                        double T) {
    std::array<double, 4> results;
    auto sqrtT = std::sqrt(T);
    auto totalvol = v * sqrtT;
    auto d1 = std::log(S / K) / totalvol + totalvol * 0.5;
    auto d2 = d1 - totalvol;
    auto pdfd1 = pdf_n(d1);
    results[0] = S * cdf_n(d1) - K * cdf_n(d2);      // price
    results[1] = S * pdfd1 * sqrtT;                  // vega
    results[2] = -pdfd1 * d2 / v;                    // vanna
    results[3] = S * pdfd1 * d1 * d2 * T / totalvol; // volga
    return results;
}

} // namespace

TEST(Paper2024Section4, base) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    std::array<double, 4> results_average{0};
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < iters; i++) {
        double S = stock_distr(generator);
        double K = stock_distr(generator);
        double v = vol_distr(generator);
        double T = time_distr(generator);

        double result = call_price(S, K, v, T);

        results_average[0] += result;
        results_average[1] += result;
        results_average[2] += result;
        results_average[3] += result;
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "base time: " << time << std::endl;

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section4, adhoc_v1) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    std::array<double, 4> results_average{0};
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    auto [S, K, v, T] = Init<4>();
    auto y = call_price(S, K, v, T);

    CalcTree t(y, S, K, v, T);

    auto dv = d<2>(v);
    auto dSv = d(S) * d(v);
    auto dvv = d(v) * d(v);

    auto dr2 = d<2>(y);
    auto bp = Tape2(dv, dSv, dvv, dr2);

    time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < iters; i++) {
        t.set(S) = stock_distr(generator);
        t.set(K) = stock_distr(generator);
        t.set(v) = vol_distr(generator);
        t.set(T) = time_distr(generator);
        t.evaluate();

        bp.reset_der();
        bp.set(dr2) = 1.;
        bp.backpropagate(t);

        results_average[0] += t.val(y);
        results_average[1] += bp.get_d(dv);
        results_average[2] += bp.get_d(dSv);
        results_average[3] += bp.get_d(dvv);
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "base time: " << time << std::endl;

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section4, adhoc_v2) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    std::array<double, 4> results_average{0};
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    auto [S, K, v, T] = Init<4>();
    auto y = call_price(S, K, v, T);
    auto vv = vega(S, K, v, T);

    CalcTree t(y, vv, S, K, v, T);

    auto dS = d(S);
    auto dv = d(v);

    auto dvega = d(vv);
    auto bp = Tape2(dvega, dS, dv);

    time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < iters; i++) {
        t.set(S) = stock_distr(generator);
        t.set(K) = stock_distr(generator);
        t.set(v) = vol_distr(generator);
        t.set(T) = time_distr(generator);
        t.evaluate();

        bp.reset_der();
        bp.set(dvega) = 1.;
        bp.backpropagate(t);

        results_average[0] += t.val(y);
        results_average[1] += t.val(vv);
        results_average[2] += bp.get_d(dS);
        results_average[3] += bp.get_d(dv);
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "base time: " << time << std::endl;

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section4, handwritten_v1) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    std::array<double, 4> results_average{0};
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < iters; i++) {
        double S = stock_distr(generator);
        double K = stock_distr(generator);
        double v = vol_distr(generator);
        double T = time_distr(generator);

        double result = call_price(S, K, v, T);
        double vega_val = vega(S, K, v, T);
        double vanna_val = vanna(S, K, v, T);
        double volga_val = volga(S, K, v, T);

        results_average[0] += result;
        results_average[1] += vega_val;
        results_average[2] += vanna_val;
        results_average[3] += volga_val;
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "base time: " << time << std::endl;

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section4, handwritten_v2) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    std::array<double, 4> results_average{0};
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < iters; i++) {
        double S = stock_distr(generator);
        double K = stock_distr(generator);
        double v = vol_distr(generator);
        double T = time_distr(generator);

        auto results = call_price_vega_vanna_volga(S, K, v, T);

        results_average[0] += results[0];
        results_average[1] += results[1];
        results_average[2] += results[2];
        results_average[3] += results[3];
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "base time: " << time << std::endl;

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

} // namespace adhoc3
