#define CODI_ChunkSize 64 // this speeds things up on runtime
#include <codi.hpp>

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

    using type = codi::RealReverse;
    using Tape = typename type::Tape;
    Tape &tape = type::getTape();

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        type S = stock_distr(generator);
        type K = stock_distr(generator);
        type v = vol_distr(generator);
        type T = time_distr(generator);

        tape.setActive();
        tape.registerInput(S);
        tape.registerInput(K);
        tape.registerInput(v);
        tape.registerInput(T);
        type y = call_price(S, K, v, T);

        tape.registerOutput(y);
        tape.setPassive();
        y.setGradient(1.0);
        tape.evaluate();

        // average values in a single Taylor expansion
        results_average[0] += y.value();
        results_average[1] += S.getGradient();
        results_average[2] += K.getGradient();
        results_average[3] += v.getGradient();
        results_average[4] += T.getGradient();

        tape.reset();
    }

    tape.reset();
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "CoDiPack order 1 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
