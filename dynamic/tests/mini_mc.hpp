#ifndef ADHOC_TEST_MINI_MC_HPP
#define ADHOC_TEST_MINI_MC_HPP

#include <tape.hpp>

#include <cmath>
#include <random>
#include <type_traits>

using adhoc_mode = adhoc::opcode<double>;
using adhoc_t = adhoc_mode::type;

// Templated compute_result with if constexpr for tape handling
template<class T>
auto
compute_result(T x1, T x2, T x3, std::size_t num_paths) -> double
{
    std::mt19937 generator(42); // same seed for reproducibility
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    // Initial calculations (before the loop)
    T y_init = x1 * x2 + cos(x3);

    // For AD types, get the tape and save position after initial calcs
    double one_over_paths = 1.0 / static_cast<double>(num_paths);

    double result = 0.0;
    for (std::size_t path = 0; path < num_paths; ++path) {
        double z1 = normal_dist(generator);
        double z2 = normal_dist(generator);

        T y_path = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);

        // add some path dependent logic to make the tapes different between
        // paths
        if (z1 > 0.) {
            y_path += x3 * z1 * z2;
        }

        result += y_path * one_over_paths;
    }

    return result;
}

template<class T>
auto
compute_result_branch(T x1, T x2, T x3, std::size_t num_paths)
{
    std::mt19937 generator(42); // same seed for reproducibility
    std::uniform_real_distribution<double> normal_dist(0.0, 1.0);

    // Initial calculations (before the loop)
    T y_init = x1 * x2;

    // For AD types, get the tape and save position after initial calcs
    [[maybe_unused]] adhoc_mode::tape_t::position_t pos2{};
    adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
    auto& tape = *tapeptr;

    double one_over_paths = 1.0 / static_cast<double>(num_paths);
    if constexpr (std::is_same_v<T, adhoc_t>) {
        pos2 = tape.get_position();
        tape.set_checkpoint();
    }

    T result = 0.0;
    for (std::size_t path = 0; path < num_paths; ++path) {
        double z1 = normal_dist(generator);
        double z2 = normal_dist(generator);

        T y_path = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);

        // add some path dependent logic to make the tapes different between
        // paths
        if (path % 2 == 0) {
            y_path += x3 * z1 * z2;
        }

        T res = y_path * one_over_paths;

        result += res;

        if constexpr (std::is_same_v<T, adhoc_t>) {
            if (tape.get_method() == adhoc::Method::Fwd) {
                tape.register_output_variable(res);
                tape.backpropagate_and_reset_to(pos2);
            }
            else {
                tape.register_output_variable(result);
                tape.set_derivative(result, 1.0);
                result.reset_id();
                tape.backpropagate_and_reset_to(pos2);
            }
        }
    }

    return result;
}

#endif // ADHOC_BUFFER_T_HPP
