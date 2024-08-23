#include <calc_tree.hpp>
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "../call_price.hpp"
#include "digital1.hpp"
#include "digital2.hpp"
#include "digital3.hpp"

#include <gtest/gtest.h>

// to make nag/dco tests work, the headers have to be copied in the folder
// ../packages/dco/. To obtain a licence contact NAG group at
// https://nag.com/contact-us/
#define NAGDCO_IN_PACKAGES_FOLDER 1

// to make CoDiPack tests work, the headers have to be copied in the folder
// ../packages/codi. CoDiPack is open source and can be obtained from
// https://github.com/SciCompKL/CoDiPack
#define CODIPACK_IN_PACKAGES_FOLDER 1

#if NAGDCO_IN_PACKAGES_FOLDER
#define DCO_AUTO_SUPPORT
#define DCO_DISABLE_AVX2_WARNING
#include "../packages/dco/include/dco.hpp"
#endif

#if CODIPACK_IN_PACKAGES_FOLDER
#include "../packages/codi/codi.hpp"
#endif

#include <chrono>
#include <random>

namespace adhoc3 {

TEST(Paper2024Section5, Primal) {
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

    std::array<double, 1> results_average{0};

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        double S = stock_distr(generator);
        double K = stock_distr(generator);
        double v = vol_distr(generator);
        double T = time_distr(generator);
        auto results = digital_price(S, K, v, T);

        // why do we do this? to make sure that compiler does not remove the
        // calculations (it might happen if they are not used)
        results_average[0] += results;
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "primal time: " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section5, AdhocOrder1) {
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

    std::array<double, 5> results_average{0};
    std::array<double, 5> results{0};

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        double S_val = stock_distr(generator);
        double K_val = stock_distr(generator);
        double v_val = vol_distr(generator);
        double T_val = time_distr(generator);
        digital_order1(S_val, K_val, v_val, T_val, results);

        // why do we do this? to make sure that compiler does not remove the
        // calculations (it might happen if they are not used)
        results_average[0] += results[0];
        results_average[1] += results[1];
        results_average[2] += results[2];
        results_average[3] += results[3];
        results_average[4] += results[4];
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "adhoc order 1 time: " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section5, AdhocOrder2) {
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

    std::array<double, 15> results_average{0};
    std::array<double, 15> results{0};

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        double S_val = stock_distr(generator);
        double K_val = stock_distr(generator);
        double v_val = vol_distr(generator);
        double T_val = time_distr(generator);
        digital_order2(S_val, K_val, v_val, T_val, results);

        // why do we do this? to make sure that compiler does not remove the
        // calculations (it might happen if they are not used)
        results_average[0] += results[0];
        results_average[1] += results[1];
        results_average[2] += results[2];
        results_average[3] += results[3];
        results_average[4] += results[4];
        results_average[5] += results[5];
        results_average[6] += results[6];
        results_average[7] += results[7];
        results_average[8] += results[8];
        results_average[9] += results[9];
        results_average[10] += results[10];
        results_average[11] += results[11];
        results_average[12] += results[12];
        results_average[13] += results[13];
        results_average[14] += results[14];
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "adhoc order 2 time: " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section5, AdhocOrder3) {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using namespace std::chrono;
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;
    double time = 0.;

    std::size_t iters = 1000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 35> results_average{0};
    std::array<double, 35> results{0};

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        double S_val = stock_distr(generator);
        double K_val = stock_distr(generator);
        double v_val = vol_distr(generator);
        double T_val = time_distr(generator);
        digital_order3(S_val, K_val, v_val, T_val, results);

        // why do we do this? to make sure that compiler does not remove the
        // calculations (it might happen if they are not used)
        results_average[0] += results[0];
        results_average[1] += results[1];
        results_average[2] += results[2];
        results_average[3] += results[3];
        results_average[4] += results[4];
        results_average[5] += results[5];
        results_average[6] += results[6];
        results_average[7] += results[7];
        results_average[8] += results[8];
        results_average[9] += results[9];
        results_average[10] += results[10];
        results_average[11] += results[11];
        results_average[12] += results[12];
        results_average[13] += results[13];
        results_average[14] += results[14];
        results_average[15] += results[15];
        results_average[16] += results[16];
        results_average[17] += results[17];
        results_average[18] += results[18];
        results_average[19] += results[19];
        results_average[20] += results[20];
        results_average[21] += results[21];
        results_average[22] += results[22];
        results_average[23] += results[23];
        results_average[24] += results[24];
        results_average[25] += results[25];
        results_average[26] += results[26];
        results_average[27] += results[27];
        results_average[28] += results[28];
        results_average[29] += results[29];
        results_average[30] += results[30];
        results_average[31] += results[31];
        results_average[32] += results[32];
        results_average[33] += results[33];
        results_average[34] += results[34];
    }

    time2 = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
               .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "adhoc order 3 time: " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

#if NAGDCO_IN_PACKAGES_FOLDER
TEST(Paper2024Section5, DcoOrder1) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average{0};

    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using mode_t = dco::ga1s<double>;
    using type = mode_t::type;
    mode_t::global_tape = mode_t::tape_t::create();

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < iters; i++) {
        type S = stock_distr(generator);
        type K = stock_distr(generator);
        type v = vol_distr(generator);
        type T = time_distr(generator);
        mode_t::global_tape->register_variable(S);
        mode_t::global_tape->register_variable(K);
        mode_t::global_tape->register_variable(v);
        mode_t::global_tape->register_variable(T);
        auto pos = mode_t::global_tape->get_position();
        type y = digital_price(S, K, v, T);
        mode_t::global_tape->register_output_variable(y);
        dco::derivative(y) = 1.0;
        mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

        // why do we do this? to make sure that compiler does not remove the
        // calculations (it might happen if they are not used)
        results_average[0] += dco::value(y);
        results_average[1] += dco::derivative(S);
        results_average[2] += dco::derivative(K);
        results_average[3] += dco::derivative(v);
        results_average[4] += dco::derivative(T);
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "nag order 1 time: " << time << std::endl;

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    mode_t::tape_t::remove(mode_t::global_tape);
}

TEST(Paper2024Section5, DcoOrder2) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 15> results_average{0};

    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using mode_t = dco::ga1s<dco::gt1s<double>::type>;
    using type = mode_t::type;
    mode_t::global_tape = mode_t::tape_t::create();

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < iters; i++) {
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
            type y = digital_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[0] += dco::value(dco::value(y));
            results_average[1] += dco::value(dco::derivative(S));
            results_average[2] += dco::value(dco::derivative(K));
            results_average[3] += dco::value(dco::derivative(v));
            results_average[4] += dco::value(dco::derivative(T));
            results_average[5] += dco::derivative(dco::derivative(S));
            results_average[9] += dco::derivative(dco::derivative(K));
            results_average[10] += dco::derivative(dco::derivative(v));
            results_average[11] += dco::derivative(dco::derivative(T));

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
            type y = digital_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[6] += dco::derivative(dco::derivative(K));
            results_average[12] += dco::derivative(dco::derivative(v));
            results_average[13] += dco::derivative(dco::derivative(T));

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
            type y = digital_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[7] += dco::derivative(dco::derivative(v));
            results_average[14] += dco::derivative(dco::derivative(T));

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
            type y = digital_price(S, K, v, T);
            mode_t::global_tape->register_output_variable(y);
            dco::derivative(y) = 1.0;
            mode_t::global_tape->interpret_adjoint_and_reset_to(pos);

            results_average[8] += dco::derivative(dco::derivative(T));

            mode_t::global_tape->reset();
            dco::derivative(dco::value(T)) = 0.0;
        }
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "nag order 2 time: " << time << std::endl;

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    mode_t::tape_t::remove(mode_t::global_tape);
}
#endif

#if CODIPACK_IN_PACKAGES_FOLDER
TEST(Paper2024Section5, CodiOrder1) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average{0};

    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using type = codi::RealReverse;
    using Tape = typename type::Tape;
    Tape &tape = type::getTape();

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < iters; i++) {
        type S = stock_distr(generator);
        type K = stock_distr(generator);
        type v = vol_distr(generator);
        type T = time_distr(generator);

        tape.setActive();
        tape.registerInput(S);
        tape.registerInput(K);
        tape.registerInput(v);
        tape.registerInput(T);
        type y = digital_price(S, K, v, T);

        tape.registerOutput(y);
        tape.setPassive();
        y.setGradient(1.0);
        tape.evaluate();

        // why do we do this? to make sure that compiler does not remove the
        // calculations (it might happen if they are not used)
        results_average[0] += y.value();
        results_average[1] += S.getGradient();
        results_average[2] += K.getGradient();
        results_average[3] += v.getGradient();
        results_average[4] += T.getGradient();
        tape.reset();
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "codipack order 1 time: " << time << std::endl;

    tape.reset();

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}

TEST(Paper2024Section5, CodiOrder2) {
    std::size_t iters = 10000;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 15> results_average{0};

    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using ts = codi::RealForwardGen<double>;
    using type = codi::RealReverseGen<ts>;
    using Tape = typename type::Tape;
    Tape &tape = type::getTape();

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < iters; i++) {
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

            type y = digital_price(S, K, v, T);
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
            S.value().gradient() = 0.0;
        }

        // 2
        {
            tape.setActive();
            tape.registerInput(S);
            tape.registerInput(K);
            tape.registerInput(v);
            tape.registerInput(T);
            K.value().gradient() = 1.0;

            type y = digital_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[6] += K.gradient().gradient();
            results_average[12] += v.gradient().gradient();
            results_average[13] += T.gradient().gradient();

            tape.reset();
            K.value().gradient() = 0.0;
        }

        // 3
        {
            tape.setActive();
            tape.registerInput(S);
            tape.registerInput(K);
            tape.registerInput(v);
            tape.registerInput(T);
            v.value().gradient() = 1.0;

            type y = digital_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[7] += v.gradient().gradient();
            results_average[14] += T.gradient().gradient();

            tape.reset();
            v.value().gradient() = 0.0;
        }

        // 4
        {
            tape.setActive();
            tape.registerInput(S);
            tape.registerInput(K);
            tape.registerInput(v);
            tape.registerInput(T);
            T.value().gradient() = 1.0;

            type y = digital_price(S, K, v, T);
            tape.registerOutput(y);
            y.setGradient(1.0);
            tape.evaluate();

            results_average[8] += T.gradient().gradient();

            tape.reset();
            T.value().gradient() = 0.0;
        }
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();
    std::cout << "iterations: " << iters << std::endl;
    std::cout << "codipack order 2 time: " << time << std::endl;

    tape.reset();

    for (std::size_t i = 0; i < results_average.size(); i++) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }
}
#endif

} // namespace adhoc3
