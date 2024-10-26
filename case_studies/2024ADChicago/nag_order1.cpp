#define DCO_AUTO_SUPPORT
#define DCO_DISABLE_AVX2_WARNING
#include "packages/dco/include/dco.hpp"

#include "black_scholes.hpp"

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

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average;
    results_average.fill(0);

    using mode_t = dco::ga1s<double>;
    using type = mode_t::type;
    mode_t::global_tape = mode_t::tape_t::create();

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        type S = stock_distr(generator);
        type K = stock_distr(generator);
        type v = vol_distr(generator);
        type T = time_distr(generator);
        mode_t::global_tape->register_variable(S);
        mode_t::global_tape->register_variable(K);
        mode_t::global_tape->register_variable(v);
        mode_t::global_tape->register_variable(T);
        auto pos = mode_t::global_tape->get_position();
        type y = call_price(S, K, v, T);
        mode_t::global_tape->register_output_variable(y);
        dco::derivative(y) = 1.0;
        mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

        // average to make sure compiler doesn't optimise calculations away
        results_average[0] += dco::value(y);
        results_average[1] += dco::derivative(S);
        results_average[2] += dco::derivative(K);
        results_average[3] += dco::derivative(v);
        results_average[4] += dco::derivative(T);
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "nag/dco order 1 time: " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    mode_t::tape_t::remove(mode_t::global_tape);
    return 0;
}
