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

    std::array<double, 15> results_average;
    results_average.fill(0);

    using ts = codi::RealForwardGen<double>;
    using type = codi::RealReverseGen<ts>;
    using Tape = typename type::Tape;
    Tape &tape = type::getTape();

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        type S = stock_distr(generator);
        type K = stock_distr(generator);
        type v = vol_distr(generator);
        type T = time_distr(generator);

        // 1
        {
            tape.setActive();
            tape.registerInput(S);
            tape.registerInput(K);
            tape.registerInput(v);
            tape.registerInput(T);
            S.value().gradient() = 1.0;

            type y = call_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[0] += y.value().value();
            results_average[1] += S.gradient().value();
            results_average[2] += K.gradient().value();
            results_average[3] += v.gradient().value();
            results_average[4] += T.gradient().value();
            results_average[5] += S.gradient().gradient();
            results_average[9] += K.gradient().gradient();
            results_average[10] += v.gradient().gradient();
            results_average[11] += T.gradient().gradient();

            tape.reset();
            // S.gradient() = 0.0;
            // K.gradient() = 0.0;
            // v.gradient() = 0.0;
            // T.gradient() = 0.0;
            S.value().gradient() = 0.0;
            // tape.reset(false);
        }

        // 2
        {
            tape.setActive();
            tape.registerInput(S);
            tape.registerInput(K);
            tape.registerInput(v);
            tape.registerInput(T);
            K.value().gradient() = 1.0;

            type y = call_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[6] += K.gradient().gradient();
            results_average[12] += v.gradient().gradient();
            results_average[13] += T.gradient().gradient();

            tape.reset();
            // S.gradient() = 0.0;
            // K.gradient() = 0.0;
            // v.gradient() = 0.0;
            // T.gradient() = 0.0;
            K.value().gradient() = 0.0;
            // tape.reset(false);
        }

        // 3
        {
            tape.setActive();
            tape.registerInput(S);
            tape.registerInput(K);
            tape.registerInput(v);
            tape.registerInput(T);
            v.value().gradient() = 1.0;

            type y = call_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[7] += v.gradient().gradient();
            results_average[14] += T.gradient().gradient();

            tape.reset();
            // S.gradient() = 0.0;
            // K.gradient() = 0.0;
            // v.gradient() = 0.0;
            // T.gradient() = 0.0;
            v.value().gradient() = 0.0;
            // tape.reset(false);
        }

        // 4
        {
            tape.setActive();
            tape.registerInput(S);
            tape.registerInput(K);
            tape.registerInput(v);
            tape.registerInput(T);
            T.value().gradient() = 1.0;

            type y = call_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[8] += T.gradient().gradient();

            tape.reset();
            // S.gradient() = 0.0;
            // K.gradient() = 0.0;
            // v.gradient() = 0.0;
            // T.gradient() = 0.0;
            T.value().gradient() = 0.0;
            // tape.reset(false);
        }
    }

    tape.reset();
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "CoDiPack order 2 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
