
#define WITHOUT_REF
#include "black_scholes.hpp"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

struct double4 {
    double S, K, v, T;
};

double4 __enzyme_autodiff(void *, ...);
int enzyme_const, enzyme_dup, enzyme_out;

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0),
        vol_distr(0.05, 0.3), time_distr(0.5, 1.5);

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average{};
    double S, K, v, T;

    auto time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        S = stock_distr(generator);
        K = stock_distr(generator);
        v = vol_distr(generator);
        T = time_distr(generator);

        auto [mu_S, mu_K, mu_v, mu_T] = __enzyme_autodiff(
            (void *)call_price<double, double, double, double>, enzyme_out, S,
            enzyme_out, K, enzyme_out, v, enzyme_out, T);

        results_average[0] += call_price(S, K, v, T);
        results_average[1] += mu_S;
        results_average[2] += mu_K;
        results_average[3] += mu_v;
        results_average[4] += mu_T;
    }

    auto time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "enzyme order 1 time (ms): " << time << std::endl;

    return 0;
}
