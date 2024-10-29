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

    std::array<double, 35> results_average{};

    auto time1 = std::chrono::high_resolution_clock::now();

    int n = 4;
    int m = 1;
    int degree = 3;
    int *multi = new int[degree];
    for (std::size_t i = 0; i < degree; ++i) {
        multi[i] = 0;
    }

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

    int dim = binomi(n + degree, degree);
    double **tensor;
    tensor = myalloc2(m, dim);

    double **S = new double *[n];

    for (std::size_t i = 0; i < n; i++) {
        S[i] = new double[n];
        for (std::size_t j = 0; j < n; j++)
            S[i][j] = (i == j) ? 1.0 : 0.0;
    }

    for (std::size_t j = 0; j < iters; ++j) {
        xp[0] = stock_distr(generator);
        xp[1] = stock_distr(generator);
        xp[2] = vol_distr(generator);
        xp[3] = time_distr(generator);

        tensor_eval(tag, m, n, degree, n, xp, tensor, S);

        // average values in a single Taylor expansion
        std::size_t counter = 0;
        for (std::size_t dim2 = 0; dim2 < (n + 1); ++dim2) {
            multi[2] = dim2;
            for (std::size_t dim1 = dim2; dim1 < (n + 1); ++dim1) {
                multi[1] = dim1;
                for (std::size_t dim0 = dim1; dim0 < (n + 1); ++dim0) {
                    multi[0] = dim0;
                    results_average[counter++] +=
                        tensor[0][tensor_address(degree, multi)];
                }
            }
        }
    }

    auto time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "ADOL-C order 3 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
