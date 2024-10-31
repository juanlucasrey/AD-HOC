#include <adolc/adolc.h>

#include "black_scholes_no_erfc.hpp"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

int factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0),
        vol_distr(0.05, 0.3), time_distr(0.5, 1.5);

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 15> results_average{};

    auto time1 = std::chrono::high_resolution_clock::now();

    int n = 4;
    int m = 1;

    double *xp = new double[n];
    double *yp = new double[m];

    auto x = new adouble[n];
    auto y = new adouble[m];

    int tag = 0;
    trace_on(tag);

    for (std::size_t i = 0; i < n; ++i) {
        x[i] <<= xp[i];
    }
    y[0] = call_price(x[0], x[1], x[2], x[3]);
    y[0] >>= yp[0];

    trace_off();

    // right multiplier of the hessian
    double *x_tangent = myalloc1(n);
    for (int i = 0; i < n; i++)
        x_tangent[i] = 0;

    // result of the forward mode pass (stores jac)
    double *y_tangent = myalloc1(m);

    // left multiplier of the hessian
    double *weights = myalloc1(m);
    weights[0] = 1.0;

    // storage for the result of hov_reverse
    double **res = myalloc2(n, 2);

    for (std::size_t j = 0; j < iters; ++j) {
        xp[0] = stock_distr(generator);
        xp[1] = stock_distr(generator);
        xp[2] = vol_distr(generator);
        xp[3] = time_distr(generator);

        std::size_t counter = 0;
        for (std::size_t i = 0; i < n; ++i) {
            x_tangent[i] = 1.0;
            fos_forward(tag, m, n, 2, xp, x_tangent, yp, y_tangent);
            hos_reverse(tag, m, n, 1, weights, res);
            results_average[i + 1] += y_tangent[0];

            // access the entries of the hessian (w/o dublication)
            for (std::size_t k = i; k < n; ++k)
                results_average[5 + counter++] += res[k][1];

            x_tangent[i] = 0.0;
        }
        results_average[0] += yp[0];
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
    std::cout << "ADOL-C fov-rev order 2 time (ms): " << time << std::endl;

    return 0;
}
