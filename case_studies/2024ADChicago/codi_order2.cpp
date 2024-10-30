#define CODI_ChunkSize 64 // this speeds things up on runtime
#define CODI_ForcedInlines 1
#include <codi.hpp>

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

    using ts = codi::RealForwardGen<double>;
    using type = codi::RealReverseGen<ts>;
    using Tape = typename type::Tape;
    Tape &tape = type::getTape();

    auto time1 = std::chrono::high_resolution_clock::now();

    type S, K, v, T;
    tape.setActive();
    tape.registerInput(S);
    tape.registerInput(K);
    tape.registerInput(v);
    tape.registerInput(T);
    auto pos = tape.getPosition();

    for (std::size_t j = 0; j < iters; ++j) {
        S.value().value() = stock_distr(generator);
        K.value().value() = stock_distr(generator);
        v.value().value() = vol_distr(generator);
        T.value().value() = time_distr(generator);

        // 1
        {
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
            results_average[6] += K.gradient().gradient();
            results_average[7] += v.gradient().gradient();
            results_average[8] += T.gradient().gradient();

            tape.resetTo(pos);
            S.gradient() = 0.0;
            K.gradient() = 0.0;
            v.gradient() = 0.0;
            T.gradient() = 0.0;
            S.value().gradient() = 0.0;
        }

        // 2
        {
            K.value().gradient() = 1.0;

            type y = call_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[9] += K.gradient().gradient();
            results_average[10] += v.gradient().gradient();
            results_average[11] += T.gradient().gradient();

            tape.resetTo(pos);
            S.gradient() = 0.0;
            K.gradient() = 0.0;
            v.gradient() = 0.0;
            T.gradient() = 0.0;
            K.value().gradient() = 0.0;
        }

        // 3
        {
            v.value().gradient() = 1.0;

            type y = call_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[12] += v.gradient().gradient();
            results_average[13] += T.gradient().gradient();

            tape.resetTo(pos);
            S.gradient() = 0.0;
            K.gradient() = 0.0;
            v.gradient() = 0.0;
            T.gradient() = 0.0;
            v.value().gradient() = 0.0;
        }

        // 4
        {
            T.value().gradient() = 1.0;

            type y = call_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[14] += T.gradient().gradient();

            tape.resetTo(pos);
            S.gradient() = 0.0;
            K.gradient() = 0.0;
            v.gradient() = 0.0;
            T.gradient() = 0.0;
            T.value().gradient() = 0.0;
        }
    }

    tape.reset();
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
    std::cout << "CoDiPack order 2 time (ms): " << time << std::endl;

    return 0;
}
