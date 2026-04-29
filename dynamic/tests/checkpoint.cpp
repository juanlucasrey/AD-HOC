#include "../../test_simple/test_simple_include.hpp"

#include <adhoc.hpp>

#include <array>
#include <cmath>
#include <random>
#include <type_traits>

using adhoc_mode = adhoc::opcode<double>;
using adhoc_t = adhoc_mode::type;

// Templated compute_result with if constexpr for tape handling
template<typename T>
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

        if constexpr (std::is_same_v<T, adhoc_t>) {
        }
        else {
            result += y_path;
        }
    }

    return result * one_over_paths;
}

// Templated compute_result with if constexpr for tape handling
template<typename T>
auto
compute_result_branch(T x1, T x2, T x3, std::size_t num_paths) -> double
{
    std::mt19937 generator(42); // same seed for reproducibility
    std::uniform_real_distribution<double> normal_dist(0.0, 1.0);

    // Initial calculations (before the loop)
    T y_init = x1 * x2;

    // For AD types, get the tape and save position after initial calcs
    [[maybe_unused]] adhoc_mode::tape_t::position_t pos2{};
    adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;

    double one_over_paths = 1.0 / static_cast<double>(num_paths);
    if constexpr (std::is_same_v<T, adhoc_t>) {
        pos2 = tape->get_position();
        tape->set_checkpoint();
    }

    double result = 0.0;
    for (std::size_t path = 0; path < num_paths; ++path) {
        double z1 = normal_dist(generator);
        double z2 = normal_dist(generator);

        T y_path = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);

        // add some path dependent logic to make the tapes different between
        // paths
        if (path % 2 == 0) {
            y_path += x3 * z1 * z2;
        }

        if constexpr (std::is_same_v<T, adhoc_t>) {
            tape->register_output_variable(y_path);
            tape->set_derivative(y_path, one_over_paths);
            tape->backpropagate_and_reset_to(pos2);
            result += y_path.get_value();
        }
        else {
            result += y_path;
        }
    }

    return result * one_over_paths;
}

void
test_checkpoint_fd_vs_ad()
{
    constexpr std::size_t num_paths = 1000;
    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double res_fd = 0.;
    double dx1_fd = 0.;
    double dx2_fd = 0.;
    double dx3_fd = 0.;

    {
        constexpr double bump = 1e-6;
        double x1 = x1_val;
        double x2 = x2_val;
        double x3 = x3_val;

        res_fd = compute_result_branch(x1, x2, x3, num_paths);

        // Central finite differences using the same templated function
        dx1_fd =
          (compute_result_branch(x1 + bump, x2, x3, num_paths) - compute_result_branch(x1 - bump, x2, x3, num_paths)) /
          (2.0 * bump);
        dx2_fd =
          (compute_result_branch(x1, x2 + bump, x3, num_paths) - compute_result_branch(x1, x2 - bump, x3, num_paths)) /
          (2.0 * bump);
        dx3_fd =
          (compute_result_branch(x1, x2, x3 + bump, num_paths) - compute_result_branch(x1, x2, x3 - bump, num_paths)) /
          (2.0 * bump);
    }

    double res_adhoc = 0.;
    double dx1_adhoc = 0.;
    double dx2_adhoc = 0.;
    double dx3_adhoc = 0.;

    {
        using adhoc_t = adhoc_t;
        // Create tape
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        res_adhoc = compute_result_branch(x1, x2, x3, num_paths);

        tape->backpropagate();

        dx1_adhoc = tape->get_derivative(x1);
        dx2_adhoc = tape->get_derivative(x2);
        dx3_adhoc = tape->get_derivative(x3);
    }

    // check that the values are close in between methods
    EXPECT_NEAR_ABS(res_adhoc, res_fd, 1e-8);
    EXPECT_NEAR_ABS(dx1_adhoc, dx1_fd, 1e-8);
    EXPECT_NEAR_ABS(dx2_adhoc, dx2_fd, 1e-8);
    EXPECT_NEAR_ABS(dx3_adhoc, dx3_fd, 1e-8);

    EXPECT_NEAR_ABS(3.421138662549827, res_adhoc, 1e-12);
    EXPECT_NEAR_ABS(1.5027371453017027, dx1_adhoc, 1e-12);
    EXPECT_NEAR_ABS(0.80747795083939378, dx2_adhoc, 1e-12);
    EXPECT_NEAR_ABS(0.11992672558600066, dx3_adhoc, 1e-12);
}

void
test_simd8_backpropagation()
{
    // Example: 8 different Monte Carlo paths computed simultaneously
    // Each lane represents a different random path
    constexpr std::size_t num_lanes = 8;

    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    // Pre-generate random numbers for 8 paths
    std::mt19937 generator(42);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    std::array<double, num_lanes> z1_vals, z2_vals;
    for (std::size_t i = 0; i < num_lanes; ++i) {
        z1_vals[i] = normal_dist(generator);
        z2_vals[i] = normal_dist(generator);
    }

    // Arrays to store results
    std::array<double, num_lanes> results;
    std::array<double, num_lanes> dx1_simd, dx2_simd, dx3_simd;

    // SIMD8 backpropagation - all 8 lanes at once
    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;

        // Set method to SIMD8
        tape->set_method(adhoc::Method::FirstOrderSimd8);
        tape->set_lanes(8);

        // Input variables
        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;
        adhoc_t x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        // Compute 8 different outputs (one per lane)
        // We'll compute each path and register as output
        std::array<adhoc_t, num_lanes> y_paths;

        adhoc_t y_init = x1 * x2 + cos(x3);

        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            double z1 = z1_vals[lane];
            double z2 = z2_vals[lane];

            y_paths[lane] = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);

            // Add path-dependent logic
            if (z1 > 0.0) {
                y_paths[lane] += x3 * z1 * z2;
            }

            results[lane] = y_paths[lane].get_value();

            // Register output and set derivative for this lane
            tape->register_output_variable(y_paths[lane]);
            tape->set_derivative(y_paths[lane], 1.0, lane);
        }

        // Single backpropagation pass for all 8 outputs
        tape->backpropagate();

        // Get derivatives for each lane
        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            dx1_simd[lane] = tape->get_derivative(x1, lane);
            dx2_simd[lane] = tape->get_derivative(x2, lane);
            dx3_simd[lane] = tape->get_derivative(x3, lane);
        }
    }

    // Verify against individual backpropagation (lane by lane)
    std::array<double, num_lanes> dx1_ref, dx2_ref, dx3_ref;

    for (std::size_t lane = 0; lane < num_lanes; ++lane) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;

        // Use standard first-order method
        tape->set_method(adhoc::Method::FirstOrderSimple);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;
        adhoc_t x3 = x3_val;

        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        double z1 = z1_vals[lane];
        double z2 = z2_vals[lane];

        adhoc_t y_init = x1 * x2 + cos(x3);
        adhoc_t y_path = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);

        if (z1 > 0.0) {
            y_path += x3 * z1 * z2;
        }

        tape->register_output_variable(y_path);
        tape->set_derivative(y_path, 1.0);
        tape->backpropagate();

        dx1_ref[lane] = tape->get_derivative(x1);
        dx2_ref[lane] = tape->get_derivative(x2);
        dx3_ref[lane] = tape->get_derivative(x3);
    }

    // Compare SIMD8 results with reference
    for (std::size_t lane = 0; lane < num_lanes; ++lane) {
        EXPECT_NEAR_ABS(dx1_simd[lane], dx1_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(dx2_simd[lane], dx2_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(dx3_simd[lane], dx3_ref[lane], 1e-10);
    }
}

void
test_simd8_monte_carlo()
{
    // More realistic Monte Carlo simulation with 8 simultaneous paths
    constexpr std::size_t num_batches = 125; // 125 batches * 8 lanes = 1000 paths
    constexpr std::size_t num_lanes = 8;
    constexpr std::size_t total_paths = num_batches * num_lanes;

    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double result_simd = 0.0;
    double dx1_simd_total = 0.0, dx2_simd_total = 0.0, dx3_simd_total = 0.0;

    // SIMD8 Monte Carlo
    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::FirstOrderSimd8);
        tape->set_lanes(8);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;
        adhoc_t x3 = x3_val;

        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        auto pos = tape->get_position();

        std::mt19937 generator(42);
        std::normal_distribution<double> normal_dist(0.0, 1.0);

        double one_over_paths = 1.0 / static_cast<double>(total_paths);

        for (std::size_t batch = 0; batch < num_batches; ++batch) {
            adhoc_t y_init = x1 * x2 + cos(x3);

            // Process 8 paths simultaneously
            std::array<adhoc_t, num_lanes> y_paths;

            for (std::size_t lane = 0; lane < num_lanes; ++lane) {
                double z1 = normal_dist(generator);
                double z2 = normal_dist(generator);

                y_paths[lane] = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);

                if (z1 > 0.0) {
                    y_paths[lane] += x3 * z1 * z2;
                }

                result_simd += y_paths[lane].get_value();

                tape->register_output_variable(y_paths[lane]);
                tape->set_derivative(y_paths[lane], one_over_paths, lane);
            }

            // Single backprop for all 8 paths
            tape->backpropagate_and_reset_to(pos);
        }

        // Get accumulated derivatives
        dx1_simd_total = tape->get_derivative(x1);
        dx2_simd_total = tape->get_derivative(x2);
        dx3_simd_total = tape->get_derivative(x3);

        // Also sum derivatives from all lanes
        for (std::size_t lane = 1; lane < num_lanes; ++lane) {
            dx1_simd_total += tape->get_derivative(x1, lane);
            dx2_simd_total += tape->get_derivative(x2, lane);
            dx3_simd_total += tape->get_derivative(x3, lane);
        }

        result_simd *= one_over_paths;
    }

    // Reference: standard sequential backpropagation
    double result_ref = 0.0;
    double dx1_ref_total = 0.0, dx2_ref_total = 0.0, dx3_ref_total = 0.0;

    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::FirstOrderSimple);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;
        adhoc_t x3 = x3_val;

        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        auto pos = tape->get_position();

        std::mt19937 generator(42);
        std::normal_distribution<double> normal_dist(0.0, 1.0);

        double one_over_paths = 1.0 / static_cast<double>(total_paths);

        for (std::size_t path = 0; path < total_paths; ++path) {
            double z1 = normal_dist(generator);
            double z2 = normal_dist(generator);

            adhoc_t y_init = x1 * x2 + cos(x3);
            adhoc_t y_path = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);

            if (z1 > 0.0) {
                y_path += x3 * z1 * z2;
            }

            result_ref += y_path.get_value();

            tape->register_output_variable(y_path);
            tape->set_derivative(y_path, one_over_paths);
            tape->backpropagate_and_reset_to(pos);
        }

        dx1_ref_total = tape->get_derivative(x1);
        dx2_ref_total = tape->get_derivative(x2);
        dx3_ref_total = tape->get_derivative(x3);

        result_ref *= one_over_paths;
    }

    EXPECT_NEAR_ABS(result_simd, result_ref, 1e-8);
    EXPECT_NEAR_ABS(dx1_simd_total, dx1_ref_total, 1e-8);
    EXPECT_NEAR_ABS(dx2_simd_total, dx2_ref_total, 1e-8);
    EXPECT_NEAR_ABS(dx3_simd_total, dx3_ref_total, 1e-8);
}

void
test_simd8_second_order_backpropagation()
{
    // Test SIMD8 second-order backpropagation
    // Each lane represents a different path with the same computation graph
    constexpr std::size_t num_lanes = 8;

    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    // Pre-generate random numbers for 8 paths
    std::mt19937 generator(42);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    std::array<double, num_lanes> z1_vals, z2_vals;
    for (std::size_t i = 0; i < num_lanes; ++i) {
        z1_vals[i] = normal_dist(generator);
        z2_vals[i] = normal_dist(generator);
    }

    // Arrays to store results
    std::array<double, num_lanes> results;
    std::array<double, num_lanes> dx1_simd, dx2_simd, dx3_simd;
    std::array<double, num_lanes> d2x1x1_simd, d2x1x2_simd, d2x1x3_simd;
    std::array<double, num_lanes> d2x2x2_simd, d2x2x3_simd, d2x3x3_simd;

    // SIMD8 second-order backpropagation - all 8 lanes at once
    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;

        // Set method to SecondOrderSimd8
        tape->set_method(adhoc::Method::SecondOrderSimd8_ankerl);
        tape->set_lanes(8);

        // Input variables
        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;
        adhoc_t x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        // Compute a simple expression: y = x1 * x2 + x3^2
        adhoc_t y = x1 * x2 + x3 * x3;

        results[0] = y.get_value();

        // Register output and set derivative for each lane
        tape->register_output_variable(y);
        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            tape->set_derivative(y, 1.0, lane);
        }

        // Single backpropagation pass for all 8 lanes
        tape->backpropagate();

        // Get first-order derivatives for each lane
        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            dx1_simd[lane] = tape->get_derivative(x1, lane);
            dx2_simd[lane] = tape->get_derivative(x2, lane);
            dx3_simd[lane] = tape->get_derivative(x3, lane);

            // Get second-order derivatives
            d2x1x1_simd[lane] = tape->get_derivative(x1, x1, lane);
            d2x1x2_simd[lane] = tape->get_derivative(x1, x2, lane);
            d2x1x3_simd[lane] = tape->get_derivative(x1, x3, lane);
            d2x2x2_simd[lane] = tape->get_derivative(x2, x2, lane);
            d2x2x3_simd[lane] = tape->get_derivative(x2, x3, lane);
            d2x3x3_simd[lane] = tape->get_derivative(x3, x3, lane);
        }
    }

    // Verify against reference (SecondOrderSimple)
    std::array<double, num_lanes> dx1_ref, dx2_ref, dx3_ref;
    std::array<double, num_lanes> d2x1x1_ref, d2x1x2_ref, d2x1x3_ref;
    std::array<double, num_lanes> d2x2x2_ref, d2x2x3_ref, d2x3x3_ref;

    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;

        // Use standard second-order method
        tape->set_method(adhoc::Method::SecondOrderSimple);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;
        adhoc_t x3 = x3_val;

        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        adhoc_t y = x1 * x2 + x3 * x3;

        tape->register_output_variable(y);
        tape->set_derivative(y, 1.0);
        tape->backpropagate();

        // Reference values (same for all lanes since computation is identical)
        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            dx1_ref[lane] = tape->get_derivative(x1);
            dx2_ref[lane] = tape->get_derivative(x2);
            dx3_ref[lane] = tape->get_derivative(x3);

            d2x1x1_ref[lane] = tape->get_derivative(x1, x1);
            d2x1x2_ref[lane] = tape->get_derivative(x1, x2);
            d2x1x3_ref[lane] = tape->get_derivative(x1, x3);
            d2x2x2_ref[lane] = tape->get_derivative(x2, x2);
            d2x2x3_ref[lane] = tape->get_derivative(x2, x3);
            d2x3x3_ref[lane] = tape->get_derivative(x3, x3);
        }
    }

    // Compare SIMD8 results with reference
    for (std::size_t lane = 0; lane < num_lanes; ++lane) {
        // First-order derivatives
        EXPECT_NEAR_ABS(dx1_simd[lane], dx1_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(dx2_simd[lane], dx2_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(dx3_simd[lane], dx3_ref[lane], 1e-10);

        // Second-order derivatives
        EXPECT_NEAR_ABS(d2x1x1_simd[lane], d2x1x1_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(d2x1x2_simd[lane], d2x1x2_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(d2x1x3_simd[lane], d2x1x3_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(d2x2x2_simd[lane], d2x2x2_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(d2x2x3_simd[lane], d2x2x3_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(d2x3x3_simd[lane], d2x3x3_ref[lane], 1e-10);
    }
}

void
test_simd8_second_order_with_transcendentals()
{
    // Test SIMD8 second-order backpropagation with transcendental functions
    constexpr std::size_t num_lanes = 8;

    double x1_val = 0.5, x2_val = 1.0;

    // Arrays for SIMD8 results
    std::array<double, num_lanes> dx1_simd, dx2_simd;
    std::array<double, num_lanes> d2x1x1_simd, d2x1x2_simd, d2x2x2_simd;

    // SIMD8 second-order backpropagation
    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::SecondOrderSimd8_ankerl);
        tape->set_lanes(8);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;

        tape->register_variable(x1);
        tape->register_variable(x2);

        // y = exp(x1) * log(x2 + 1) + sqrt(x1 + x2)
        adhoc_t y = exp(x1) * log(x2 + 1.0) + sqrt(x1 + x2);

        tape->register_output_variable(y);
        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            tape->set_derivative(y, 1.0, lane);
        }

        tape->backpropagate();

        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            dx1_simd[lane] = tape->get_derivative(x1, lane);
            dx2_simd[lane] = tape->get_derivative(x2, lane);
            d2x1x1_simd[lane] = tape->get_derivative(x1, x1, lane);
            d2x1x2_simd[lane] = tape->get_derivative(x1, x2, lane);
            d2x2x2_simd[lane] = tape->get_derivative(x2, x2, lane);
        }
    }

    // Reference with SecondOrderSimple
    double dx1_ref, dx2_ref;
    double d2x1x1_ref, d2x1x2_ref, d2x2x2_ref;

    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::SecondOrderSimple);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;

        tape->register_variable(x1);
        tape->register_variable(x2);

        adhoc_t y = exp(x1) * log(x2 + 1.0) + sqrt(x1 + x2);

        tape->register_output_variable(y);
        tape->set_derivative(y, 1.0);
        tape->backpropagate();

        dx1_ref = tape->get_derivative(x1);
        dx2_ref = tape->get_derivative(x2);
        d2x1x1_ref = tape->get_derivative(x1, x1);
        d2x1x2_ref = tape->get_derivative(x1, x2);
        d2x2x2_ref = tape->get_derivative(x2, x2);
    }

    // Compare
    for (std::size_t lane = 0; lane < num_lanes; ++lane) {
        EXPECT_NEAR_ABS(dx1_simd[lane], dx1_ref, 1e-10);
        EXPECT_NEAR_ABS(dx2_simd[lane], dx2_ref, 1e-10);
        EXPECT_NEAR_ABS(d2x1x1_simd[lane], d2x1x1_ref, 1e-10);
        EXPECT_NEAR_ABS(d2x1x2_simd[lane], d2x1x2_ref, 1e-10);
        EXPECT_NEAR_ABS(d2x2x2_simd[lane], d2x2x2_ref, 1e-10);
    }

    // Verify against finite differences
    constexpr double bump = 1e-5;

    auto f = [](double x1, double x2) { return std::exp(x1) * std::log(x2 + 1.0) + std::sqrt(x1 + x2); };

    double dx1_fd = (f(x1_val + bump, x2_val) - f(x1_val - bump, x2_val)) / (2.0 * bump);
    double dx2_fd = (f(x1_val, x2_val + bump) - f(x1_val, x2_val - bump)) / (2.0 * bump);

    double d2x1x1_fd = (f(x1_val + bump, x2_val) - 2.0 * f(x1_val, x2_val) + f(x1_val - bump, x2_val)) / (bump * bump);
    double d2x2x2_fd = (f(x1_val, x2_val + bump) - 2.0 * f(x1_val, x2_val) + f(x1_val, x2_val - bump)) / (bump * bump);
    double d2x1x2_fd = (f(x1_val + bump, x2_val + bump) - f(x1_val + bump, x2_val - bump) -
                        f(x1_val - bump, x2_val + bump) + f(x1_val - bump, x2_val - bump)) /
                       (4.0 * bump * bump);

    EXPECT_NEAR_ABS(dx1_ref, dx1_fd, 1e-6);
    EXPECT_NEAR_ABS(dx2_ref, dx2_fd, 1e-6);
    EXPECT_NEAR_ABS(d2x1x1_ref, d2x1x1_fd, 1e-4);
    EXPECT_NEAR_ABS(d2x1x2_ref, d2x1x2_fd, 1e-4);
    EXPECT_NEAR_ABS(d2x2x2_ref, d2x2x2_fd, 1e-4);
}

void
test_simd8_second_order_monte_carlo()
{
    // Monte Carlo simulation with SIMD8 second-order derivatives
    constexpr std::size_t num_batches = 125; // 125 batches * 8 lanes = 1000 paths
    constexpr std::size_t num_lanes = 8;
    constexpr std::size_t total_paths = num_batches * num_lanes;

    double x1_val = 1.5, x2_val = 2.0;

    double result_simd = 0.0;
    double dx1_simd_total = 0.0, dx2_simd_total = 0.0;
    double d2x1x1_simd_total = 0.0, d2x1x2_simd_total = 0.0, d2x2x2_simd_total = 0.0;

    // SIMD8 second-order Monte Carlo
    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::SecondOrderSimd8_ankerl);
        tape->set_lanes(8);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;

        tape->register_variable(x1);
        tape->register_variable(x2);

        auto pos = tape->get_position();

        std::mt19937 generator(42);
        std::normal_distribution<double> normal_dist(0.0, 1.0);

        double one_over_paths = 1.0 / static_cast<double>(total_paths);

        for (std::size_t batch = 0; batch < num_batches; ++batch) {
            // Process 8 paths simultaneously
            std::array<adhoc_t, num_lanes> y_paths;

            for (std::size_t lane = 0; lane < num_lanes; ++lane) {
                double z = normal_dist(generator);

                // Simple payoff: (x1 * z + x2)^2
                adhoc_t temp = x1 * z + x2;
                y_paths[lane] = temp * temp;

                result_simd += y_paths[lane].get_value();

                tape->register_output_variable(y_paths[lane]);
                tape->set_derivative(y_paths[lane], one_over_paths, lane);
            }

            // Single backprop for all 8 paths
            tape->backpropagate_and_reset_to(pos);
        }

        // Get accumulated derivatives from all lanes
        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            dx1_simd_total += tape->get_derivative(x1, lane);
            dx2_simd_total += tape->get_derivative(x2, lane);
            d2x1x1_simd_total += tape->get_derivative(x1, x1, lane);
            d2x1x2_simd_total += tape->get_derivative(x1, x2, lane);
            d2x2x2_simd_total += tape->get_derivative(x2, x2, lane);
        }

        result_simd *= one_over_paths;
    }

    // Reference: standard sequential second-order backpropagation
    double result_ref = 0.0;
    double dx1_ref_total = 0.0, dx2_ref_total = 0.0;
    double d2x1x1_ref_total = 0.0, d2x1x2_ref_total = 0.0, d2x2x2_ref_total = 0.0;

    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::SecondOrderSimple);

        adhoc_t x1 = x1_val;
        adhoc_t x2 = x2_val;

        tape->register_variable(x1);
        tape->register_variable(x2);

        auto pos = tape->get_position();

        std::mt19937 generator(42);
        std::normal_distribution<double> normal_dist(0.0, 1.0);

        double one_over_paths = 1.0 / static_cast<double>(total_paths);

        for (std::size_t path = 0; path < total_paths; ++path) {
            double z = normal_dist(generator);

            adhoc_t temp = x1 * z + x2;
            adhoc_t y = temp * temp;

            result_ref += y.get_value();

            tape->register_output_variable(y);
            tape->set_derivative(y, one_over_paths);
            tape->backpropagate_and_reset_to(pos);
        }

        dx1_ref_total = tape->get_derivative(x1);
        dx2_ref_total = tape->get_derivative(x2);
        d2x1x1_ref_total = tape->get_derivative(x1, x1);
        d2x1x2_ref_total = tape->get_derivative(x1, x2);
        d2x2x2_ref_total = tape->get_derivative(x2, x2);

        result_ref *= one_over_paths;
    }
    EXPECT_NEAR_ABS(result_simd, result_ref, 1e-8);
    EXPECT_NEAR_ABS(dx1_simd_total, dx1_ref_total, 1e-8);
    EXPECT_NEAR_ABS(dx2_simd_total, dx2_ref_total, 1e-8);
    EXPECT_NEAR_ABS(d2x1x1_simd_total, d2x1x1_ref_total, 1e-8);
    EXPECT_NEAR_ABS(d2x1x2_simd_total, d2x1x2_ref_total, 1e-8);
    EXPECT_NEAR_ABS(d2x2x2_simd_total, d2x2x2_ref_total, 1e-8);
}

void
test_simd8_second_order_different_seeds_per_lane()
{
    // Test with different computations per lane (different random seeds)
    constexpr std::size_t num_lanes = 8;

    double x_val = 1.5;

    // Different multipliers per lane to simulate different paths
    std::array<double, num_lanes> multipliers = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 };

    // Arrays to store SIMD8 results
    std::array<double, num_lanes> dx_simd, d2xx_simd;

    // SIMD8 computation
    {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::SecondOrderSimd8_ankerl);
        tape->set_lanes(8);

        adhoc_t x = x_val;
        tape->register_variable(x);

        // Compute y = (x * m)^2 where m is different per lane
        // But we can only have one tape, so we compute with m=1 and scale
        // derivatives
        adhoc_t y = x * x;

        tape->register_output_variable(y);
        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            // Scale the seed derivative by multiplier^2 to account for chain
            // rule
            tape->set_derivative(y, multipliers[lane] * multipliers[lane], lane);
        }

        tape->backpropagate();

        for (std::size_t lane = 0; lane < num_lanes; ++lane) {
            dx_simd[lane] = tape->get_derivative(x, lane);
            d2xx_simd[lane] = tape->get_derivative(x, x, lane);
        }
    }

    // Reference: compute each lane separately
    std::array<double, num_lanes> dx_ref, d2xx_ref;

    for (std::size_t lane = 0; lane < num_lanes; ++lane) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::SecondOrderSimple);

        adhoc_t x = x_val;
        tape->register_variable(x);

        adhoc_t y = x * x;

        tape->register_output_variable(y);
        tape->set_derivative(y, multipliers[lane] * multipliers[lane]);
        tape->backpropagate();

        dx_ref[lane] = tape->get_derivative(x);
        d2xx_ref[lane] = tape->get_derivative(x, x);
    }

    // Compare
    for (std::size_t lane = 0; lane < num_lanes; ++lane) {
        EXPECT_NEAR_ABS(dx_simd[lane], dx_ref[lane], 1e-10);
        EXPECT_NEAR_ABS(d2xx_simd[lane], d2xx_ref[lane], 1e-10);
    }
}

void
test_checkpoint_branch_lossy()
{
    constexpr std::size_t num_paths = 1000;
    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double res_bwd = 0.;
    double dx1_bwd = 0.;
    double dx2_bwd = 0.;
    double dx3_bwd = 0.;

    {
        using adhoc_t = adhoc_t;
        // Create tape
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        res_bwd = compute_result_branch(x1, x2, x3, num_paths);

        // tape->set_derivative(res_bwd, 1.0);
        tape->backpropagate();
        // adhoc::BackPropagator2 bp(tape);
        // bp.set_derivative(res_bwd, 1.0);
        // bp.backpropagate();

        // Retrieve accumulated derivatives
        dx1_bwd = tape->get_derivative(x1);
        dx2_bwd = tape->get_derivative(x2);
        dx3_bwd = tape->get_derivative(x3);

        // Basic sanity check - derivatives should be non-zero
        EXPECT_NOT_EQUAL(dx1_bwd, 0.0);
        EXPECT_NOT_EQUAL(dx2_bwd, 0.0);
        EXPECT_NOT_EQUAL(dx3_bwd, 0.0);
    }

    double res_fd = 0.;
    double dx1_fd = 0.;
    double dx2_fd = 0.;
    double dx3_fd = 0.;

    {
        constexpr double bump = 1e-6;
        double x1 = x1_val;
        double x2 = x2_val;
        double x3 = x3_val;

        res_fd = compute_result_branch(x1, x2, x3, num_paths);

        // Central finite differences using the same templated function
        dx1_fd =
          (compute_result_branch(x1 + bump, x2, x3, num_paths) - compute_result_branch(x1 - bump, x2, x3, num_paths)) /
          (2.0 * bump);
        dx2_fd =
          (compute_result_branch(x1, x2 + bump, x3, num_paths) - compute_result_branch(x1, x2 - bump, x3, num_paths)) /
          (2.0 * bump);
        dx3_fd =
          (compute_result_branch(x1, x2, x3 + bump, num_paths) - compute_result_branch(x1, x2, x3 - bump, num_paths)) /
          (2.0 * bump);
    }

    // check that the values are close in between methods
    EXPECT_NEAR_ABS(res_bwd, res_fd, 1e-8);
    EXPECT_NEAR_ABS(dx1_bwd, dx1_fd, 1e-8);
    EXPECT_NEAR_ABS(dx2_bwd, dx2_fd, 1e-8);
    EXPECT_NEAR_ABS(dx3_bwd, dx3_fd, 1e-8);

    double res_lossy = 0.;
    double dx1_lossy = 0.;
    double dx2_lossy = 0.;
    double dx3_lossy = 0.;

    {
        using adhoc_t = adhoc_t;
        // Create tape
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(adhoc::Method::FirstOrderLossy);

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        res_lossy = compute_result_branch(x1, x2, x3, num_paths);

        // tape->set_derivative(res_bwd, 1.0);
        tape->backpropagate();
        // adhoc::BackPropagator2 bp(tape);
        // bp.set_derivative(res_bwd, 1.0);
        // bp.backpropagate();

        // Retrieve accumulated derivatives
        dx1_lossy = tape->get_derivative(x1);
        dx2_lossy = tape->get_derivative(x2);
        dx3_lossy = tape->get_derivative(x3);

        // Basic sanity check - derivatives should be non-zero
        EXPECT_NOT_EQUAL(dx1_lossy, 0.0);
        EXPECT_NOT_EQUAL(dx2_lossy, 0.0);
        EXPECT_NOT_EQUAL(dx3_lossy, 0.0);
    }

    EXPECT_NEAR_ABS(res_bwd, res_lossy, 1e-8);
    EXPECT_NEAR_ABS(dx1_bwd, dx1_lossy, 1e-8);
    EXPECT_NEAR_ABS(dx2_bwd, dx2_lossy, 1e-8);
    EXPECT_NEAR_ABS(dx3_bwd, dx3_lossy, 1e-8);
}

void
test_lossy_compressed_1in1out()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };
    double deriv_val = -199.95882307798402;
    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);
        adhoc_t x = 2.3;
        tape->register_variable(x);
        auto y = log(erfc(exp(x)));
        tape->register_output_variable(y);
        tape->set_derivative(y, 1.);
        tape->backpropagate();
        EXPECT_NEAR_ABS(deriv_val, tape->get_derivative(x), 1e-8);
    }
}

void
test_lossy_compressed_2in1out()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };
    double deriv_val1 = 518.81914462776513;
    double deriv_val2 = -269.97895210733441;
    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);
        adhoc_t x1 = 2.3;
        adhoc_t x2 = 1.1;
        tape->register_variable(x1);
        tape->register_variable(x2);

        auto y1 = log(erfc(exp(x1)));
        auto y2 = exp(erfc(log(x2)));
        auto y = erf(cos(y1 * y2));

        tape->register_output_variable(y);
        tape->set_derivative(y, 1.);
        tape->backpropagate();
        EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-8);
        EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x2), 1e-8);
    }
}

void
test_lossy_compressed_1in1out2paths()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };
    double deriv_val1 = -2.381398009286281;
    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);
        adhoc_t x1 = 2.3;
        tape->register_variable(x1);

        auto y = cos(x1) * x1;

        tape->register_output_variable(y);
        tape->set_derivative(y, 1.);
        tape->backpropagate();
        EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-13);
    }
}

void
test_lossy_compressed_1in1out2paths2()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };
    double deriv_val1 = -2.381398009286281;
    double deriv_val2 = -5.4772154213584461;
    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);
        adhoc_t x1 = 2.3;
        adhoc_t x2 = 1;
        tape->register_variable(x1);
        tape->register_variable(x2);

        auto prod = x1 * x2;
        auto y = cos(prod) * prod;

        tape->register_output_variable(y);
        tape->set_derivative(y, 1.);
        tape->backpropagate();

        EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-13);
        EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x2), 1e-13);
    }
}

void
test_lossy_compressed_2in1out2()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };
    double deriv_val1 = 42.516292658681806;
    double deriv_val2 = -22.124287930267304;
    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);
        adhoc_t x1 = 2.3;
        adhoc_t x2 = 1.1;
        tape->register_variable(x1);
        tape->register_variable(x2);

        auto y1 = log(erfc(exp(x1)));
        auto y2 = exp(erfc(log(x2)));
        auto prod = y1 * y2;
        auto y = erf(cos(prod)) * exp(prod * 0.01);

        tape->register_output_variable(y);
        tape->set_derivative(y, 1.);
        tape->backpropagate();

        EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-8);
        EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x2), 1e-8);
    }
}

void
test_lossy_compressed_1in2out()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };

    double deriv_val1 = 752.01628170139816;
    double deriv_val2 = -0.0027146390620349445;
    double deriv_val3 = 376.00678353116808;
    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        auto& tape_ref = *tape;
        tape->set_method(m);
        adhoc_t x1 = 2.3;
        tape->register_variable(x1);

        auto int1 = log(erfc(exp(x1)));
        auto y1 = log(cos(int1));
        auto y2 = erfc(exp(int1 * -0.01));

        tape->register_output_variable(y1);
        tape->register_output_variable(y2);
        tape->set_derivative(y1, 1.);
        tape->set_derivative(y2, 0.);
        tape->backpropagate();

        EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-8);

        tape->zero_adjoints();
        tape->set_derivative(y1, 0.);
        tape->set_derivative(y2, 1.);
        tape->backpropagate();

        EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x1), 1e-8);

        tape->zero_adjoints();
        tape->set_derivative(y1, 0.5);
        tape->set_derivative(y2, 0.5);
        tape->backpropagate();

        EXPECT_NEAR_ABS(deriv_val3, tape->get_derivative(x1), 1e-8);

        // if (m == adhoc::Method::FirstOrderSimple) {
        //     deriv_val = tape->get_derivative(x1);
        // }
        // else {
        //     EXPECT_NEAR_ABS(deriv_val, tape->get_derivative(x1), 1e-13);
        // }
    }
}

void
test_lossy_compressed_complex1_pre()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };

    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double result_val = 0.;
    double deriv_val1 = 0.;
    double deriv_val2 = 0.;
    double deriv_val3 = 0.;

    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        double z1 = 0.8;
        double z2 = 0.2;

        auto y_init = x1 * x2;
        auto y_path = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);
        auto y = y_path + x3 * z1 * z2;
        tape->register_output_variable(y);
        tape->set_derivative(y, 1.0);
        tape->backpropagate();

        if (m == adhoc::Method::FirstOrderSimple) {
            result_val = y.get_value();
            deriv_val1 = tape->get_derivative(x1);
            deriv_val2 = tape->get_derivative(x2);
            deriv_val3 = tape->get_derivative(x3);
        }
        else {
            EXPECT_NEAR_ABS(result_val, y.get_value(), 1e-13);
            EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-13);
            EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x2), 1e-13);
            EXPECT_NEAR_ABS(deriv_val3, tape->get_derivative(x3), 1e-13);
        }
    }
}

void
test_lossy_compressed_complex1_pre2()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };

    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double result_val = 0.;
    double deriv_val1 = 0.;
    double deriv_val2 = 0.;
    double deriv_val3 = 0.;

    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        // Create tape
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        auto y_init = x1 * x2;

        double one_over_paths = 1.0;

        auto pos2 = tape->get_position();
        tape->set_checkpoint();

        double result = 0.0;
        double z1 = 0.8;
        double z2 = 0.2;

        auto y_path = y_init * 0.8;

        tape->register_output_variable(y_path);
        tape->set_derivative(y_path, one_over_paths);
        tape->backpropagate_and_reset_to(pos2);
        result += y_path.get_value();

        double res_lossy = result * one_over_paths;

        tape->backpropagate();

        if (m == adhoc::Method::FirstOrderSimple) {
            result_val = res_lossy;
            deriv_val1 = tape->get_derivative(x1);
            deriv_val2 = tape->get_derivative(x2);
            deriv_val3 = tape->get_derivative(x3);
        }
        else {
            EXPECT_NEAR_ABS(result_val, res_lossy, 1e-13);
            EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-13);
            EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x2), 1e-13);
            EXPECT_NEAR_ABS(deriv_val3, tape->get_derivative(x3), 1e-13);
        }

        // EXPECT_NEAR_ABS(result_val, res_lossy, 1e-13);
        // EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-13);
        // EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x2), 1e-13);
        // EXPECT_NEAR_ABS(deriv_val3, tape->get_derivative(x3), 1e-13);
    }
}

void
test_lossy_compressed_complex1()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderSimple,
                                           adhoc::Method::FirstOrderLossy,
                                           adhoc::Method::FirstOrderLossyCompressed };

    constexpr std::size_t num_paths = 1000;
    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double result_val = 3.421138662549827;
    double deriv_val1 = 1.5027371453017027;
    double deriv_val2 = 0.80747795083939388;
    double deriv_val3 = 0.11992672558600065;

    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        // Create tape
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tape;
        tape->set_method(m);

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape->register_variable(x1);
        tape->register_variable(x2);
        tape->register_variable(x3);

        double res_lossy = compute_result_branch(x1, x2, x3, num_paths);

        tape->backpropagate();

        EXPECT_NEAR_ABS(result_val, res_lossy, 1e-13);
        EXPECT_NEAR_ABS(deriv_val1, tape->get_derivative(x1), 1e-13);
        EXPECT_NEAR_ABS(deriv_val2, tape->get_derivative(x2), 1e-13);
        EXPECT_NEAR_ABS(deriv_val3, tape->get_derivative(x3), 1e-13);
    }
}

void
test_lossy_compressed_complex2()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderLossy, adhoc::Method::FirstOrderLossyCompressed };

    double df = 0;
    double f = 0;

    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;

    for (auto m : methods) {
        D inputD = 0.5;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->set_method(m);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);

        auto const scEnd = inputD * 1.2;
        auto outputD = scEnd * (scEnd * 1.1);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();

        if (m == adhoc::Method::FirstOrderLossy) {
            df = tape.get_derivative(inputD);
            f = adhoc::passive_value(outputD);
        }
        else {
            EXPECT_NEAR_ABS(df, tape.get_derivative(inputD), 1e-13);
            EXPECT_NEAR_ABS(f, adhoc::passive_value(outputD), 1e-13);
        }
    }
}

void
test_lossy_compressed_complex3()
{
    std::vector<adhoc::Method> methods = { adhoc::Method::FirstOrderLossy, adhoc::Method::FirstOrderLossyCompressed };

    double df = 0;
    double f = 0;

    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;

    for (auto m : methods) {
        D inputD = 0.5;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->set_method(m);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);

        auto const scEnd = inputD * 1.2;
        auto intermediary = scEnd * (scEnd * 1.1);
        auto outputD = exp(cos(scEnd)) * log(intermediary);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();

        if (m == adhoc::Method::FirstOrderLossy) {
            df = tape.get_derivative(inputD);
            f = adhoc::passive_value(outputD);
        }
        else {
            EXPECT_NEAR_ABS(df, tape.get_derivative(inputD), 1e-13);
            EXPECT_NEAR_ABS(f, adhoc::passive_value(outputD), 1e-13);
        }
    }
}

auto
main() -> int
{
    test_checkpoint_fd_vs_ad();
    test_checkpoint_branch_lossy();
    test_simd8_backpropagation();
    test_simd8_monte_carlo();
    test_simd8_second_order_backpropagation();
    test_simd8_second_order_with_transcendentals();
    test_simd8_second_order_monte_carlo();
    test_simd8_second_order_different_seeds_per_lane();

    test_lossy_compressed_complex2();
    test_lossy_compressed_complex3();
    test_lossy_compressed_1in1out();
    test_lossy_compressed_2in1out();
    test_lossy_compressed_1in1out2paths();
    test_lossy_compressed_1in1out2paths2();
    test_lossy_compressed_2in1out2();
    test_lossy_compressed_1in2out();
    test_lossy_compressed_complex1_pre();
    test_lossy_compressed_complex1_pre2();
    test_lossy_compressed_complex1();

    TEST_END;
}
