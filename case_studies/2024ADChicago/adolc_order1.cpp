#include <adolc/adolc.h>

#include "black_scholes_no_erfc.hpp"

#include <chrono>
#include <iostream>
#include <random>

int factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using namespace std::chrono;
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average;
    results_average.fill(0);

    time1 = std::chrono::high_resolution_clock::now();

    // Define the order of the highest derivative
    // int degree = 1;
    int n = 4;
    int m = 1;
    auto x = new adouble[n];
    auto y = new adouble[m];

    double *xp = new double[n];
    double *yp = new double[m];

    double *u = new double[m];
    double *z = new double[n];
    u[0] = 1.;

    int tag = 0;
    trace_on(tag);

    for (std::size_t i = 0; i < n; ++i) {
        x[i] <<= xp[i];
    }
    y[0] = call_price(x[0], x[1], x[2], x[3]);
    y[0] >>= yp[0];

    trace_off();

    for (std::size_t j = 0; j < iters; ++j) {
        xp[0] = stock_distr(generator);
        xp[1] = stock_distr(generator);
        xp[2] = vol_distr(generator);
        xp[3] = time_distr(generator);

        zos_forward(tag, m, n, 1, xp, yp);
        fos_reverse(tag, m, n, u, z);

        // average to make sure compiler doesn't optimise calculations away
        results_average[0] += yp[0];
        results_average[1] += z[0];
        results_average[2] += z[1];
        results_average[3] += z[2];
        results_average[4] += z[3];
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "adhoc order 1 time: " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
