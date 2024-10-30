#define DCO_AUTO_SUPPORT
#define DCO_DISABLE_AVX2_WARNING
#include <dco.hpp>

#include "black_scholes.hpp"

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

    std::array<double, 15> results_average{};

    using mode_t = dco::ga1s<dco::gt1s<double>::type>;
    using type = mode_t::type;
    mode_t::global_tape = mode_t::tape_t::create();

    auto time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        type S = stock_distr(generator);
        type K = stock_distr(generator);
        type v = vol_distr(generator);
        type T = time_distr(generator);

        // 1
        {
            mode_t::global_tape->register_variable(S);
            mode_t::global_tape->register_variable(K);
            mode_t::global_tape->register_variable(v);
            mode_t::global_tape->register_variable(T);
            auto pos = mode_t::global_tape->get_position();

            dco::derivative(dco::value(S)) = 1.0;
            type y = call_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[0] += dco::value(dco::value(y));
            results_average[1] += dco::value(dco::derivative(S));
            results_average[2] += dco::value(dco::derivative(K));
            results_average[3] += dco::value(dco::derivative(v));
            results_average[4] += dco::value(dco::derivative(T));
            results_average[5] += dco::derivative(dco::derivative(S));
            results_average[6] += dco::derivative(dco::derivative(K));
            results_average[7] += dco::derivative(dco::derivative(v));
            results_average[8] += dco::derivative(dco::derivative(T));

            mode_t::global_tape->reset();
            dco::derivative(dco::value(S)) = 0.0;
        }

        // 2
        {
            mode_t::global_tape->register_variable(S);
            mode_t::global_tape->register_variable(K);
            mode_t::global_tape->register_variable(v);
            mode_t::global_tape->register_variable(T);
            auto pos = mode_t::global_tape->get_position();

            dco::derivative(dco::value(K)) = 1.0;
            type y = call_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[9] += dco::derivative(dco::derivative(K));
            results_average[10] += dco::derivative(dco::derivative(v));
            results_average[11] += dco::derivative(dco::derivative(T));

            mode_t::global_tape->reset();
            dco::derivative(dco::value(K)) = 0.0;
        }

        // 3
        {
            mode_t::global_tape->register_variable(S);
            mode_t::global_tape->register_variable(K);
            mode_t::global_tape->register_variable(v);
            mode_t::global_tape->register_variable(T);
            auto pos = mode_t::global_tape->get_position();

            dco::derivative(dco::value(v)) = 1.0;
            type y = call_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[12] += dco::derivative(dco::derivative(v));
            results_average[13] += dco::derivative(dco::derivative(T));

            mode_t::global_tape->reset();
            dco::derivative(dco::value(v)) = 0.0;
        }

        // 4
        {
            mode_t::global_tape->register_variable(S);
            mode_t::global_tape->register_variable(K);
            mode_t::global_tape->register_variable(v);
            mode_t::global_tape->register_variable(T);
            auto pos = mode_t::global_tape->get_position();

            dco::derivative(dco::value(T)) = 1.0;
            type y = call_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[14] += dco::derivative(dco::derivative(T));

            mode_t::global_tape->reset();
            dco::derivative(dco::value(T)) = 0.0;
        }
    }

    mode_t::tape_t::remove(mode_t::global_tape);
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
    std::cout << "dco/c++ order 2 time (ms): " << time << std::endl;

    return 0;
}
