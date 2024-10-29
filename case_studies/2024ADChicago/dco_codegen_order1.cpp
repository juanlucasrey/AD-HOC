#include <array>
#include <chrono>
#include <iostream>
#include <random>

#include "dco_codegen/dco_codegen_adjoint.hpp"

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0),
        vol_distr(0.05, 0.3), time_distr(0.5, 1.5);

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average{};

    std::vector<double> inputs(4), dinputs(4, 0.0), outputs(1), doutputs{1.0},
        p;
    dco::codegen_buffer_t buffer;

    auto time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t j = 0; j < iters; ++j) {
        inputs[0] = stock_distr(generator);
        inputs[1] = stock_distr(generator);
        inputs[2] = vol_distr(generator);
        inputs[3] = time_distr(generator);

        adjoint(inputs, dinputs, p, outputs, doutputs, buffer);

        results_average[0] += outputs[0];
    }

    results_average[1] = dinputs[0];
    results_average[2] = dinputs[1];
    results_average[3] = dinputs[2];
    results_average[4] = dinputs[3];

    auto time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "dco Codegen order 1 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
