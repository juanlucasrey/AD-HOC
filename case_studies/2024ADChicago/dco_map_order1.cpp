#include <dco_map.hpp>

#include "black_scholes.hpp"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

namespace {

// template <class D>
// inline void call_price_wrap(const D &S, const D &K, const D &v, const D &T,
//                             D &y) {
//     y = call_price(S, K, v, T);
// }

template <class D>
inline void call_price_wrap(const D &S, const D &K, const D &v, const D &T,
                            D &y) {
    using adhoc::constants::CD;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<0.5>();
    auto d2 = d1 - totalvol;
    y = S * cdf_n(d1) - K * cdf_n(d2);
}

} // namespace

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0),
        vol_distr(0.05, 0.3), time_distr(0.5, 1.5);

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 5> results_average{};

    using mode_t = dco_map::ga1s<double>;
    using type = mode_t::type;

    auto time1 = std::chrono::high_resolution_clock::now();

    type S, K, v, T;
    for (std::size_t j = 0; j < iters; ++j) {
        dco_map::value(S) = stock_distr(generator);
        dco_map::value(K) = stock_distr(generator);
        dco_map::value(v) = vol_distr(generator);
        dco_map::value(T) = time_distr(generator);
        type y;
        dco_map::derivative(y) = 1.0;

        call_price_wrap(S, K, v, T, y);

        results_average[0] += dco_map::value(y);
    }

    results_average[1] = dco_map::derivative(S);
    results_average[2] = dco_map::derivative(K);
    results_average[3] = dco_map::derivative(v);
    results_average[4] = dco_map::derivative(T);
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
    std::cout << "dco/map order 1 time (ms): " << time << std::endl;

    return 0;
}
