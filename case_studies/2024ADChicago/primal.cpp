#include "black_scholes.hpp"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using namespace std::chrono;
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    std::size_t iters = 1000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 1> results_average;
    results_average.fill(0);

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        double S = stock_distr(generator);
        double K = stock_distr(generator);
        double v = vol_distr(generator);
        double T = time_distr(generator);

        // why do we do this? to make sure that compiler does not remove the
        // calculations (it might happen if they are not used)
        results_average[0] += call_price(S, K, v, T);
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "primal value time: " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
