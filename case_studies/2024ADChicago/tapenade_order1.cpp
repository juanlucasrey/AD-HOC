#include "black_scholes_c.h"
#include "tapenade/c_version_bs_b.h"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0),
        vol_distr(0.05, 0.3), time_distr(0.5, 1.5);

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average{};

    double S0, K0, v0, T0, Sb = 0, Kb = 0, vb = 0, Tb = 0;
    double seed = 1.0;

    auto time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        S0 = stock_distr(generator);
        K0 = stock_distr(generator);
        v0 = vol_distr(generator);
        T0 = time_distr(generator);
        results_average[0] += call_price(S0, K0, v0, T0);
        call_price_b(S0, &Sb, K0, &Kb, v0, &vb, T0, &Tb, seed);

        results_average[1] += Sb;
        results_average[2] += Kb;
        results_average[3] += vb;
        results_average[4] += Tb;

        Sb = 0;
        Kb = 0;
        vb = 0;
        Tb = 0;
    }

    auto time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "Tapenade order 1 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
