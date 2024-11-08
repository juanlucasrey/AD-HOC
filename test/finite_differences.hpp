#ifndef FINITE_DIFFERENCES_HPP
#define FINITE_DIFFERENCES_HPP

#include <array>
#include <cmath>
#include <functional>
#include <numeric>

template <class D>
std::array<D, 6> finite_differences(D x, D epsilon, std::function<D(D)> func) {
    std::array<D, 11> fvals{
        func(x - 5. * epsilon), func(x - 4. * epsilon), func(x - 3. * epsilon),
        func(x - 2. * epsilon), func(x - 1. * epsilon), func(x),
        func(x + 1. * epsilon), func(x + 2. * epsilon), func(x + 3. * epsilon),
        func(x + 4. * epsilon), func(x + 5. * epsilon)};

    constexpr std::array<D, 11> weights1 = {
        0.,      1. / 280., -4. / 105., 1. / 5.,   -4. / 5., 0.,
        4. / 5., -1. / 5.,  4. / 105.,  -1. / 280, 0.};

    constexpr std::array<D, 11> weights2 = {
        0.,      -1. / 560., 8. / 315., -1. / 5.,   8. / 5., -205. / 72.,
        8. / 5., -1. / 5.,   8. / 315., -1. / 560., 0.};

    constexpr std::array<D, 11> weights3 = {
        0.,         -7. / 240.,  3. / 10.,  -169. / 120., 61. / 30., 0.,
        -61. / 30., 169. / 120., -3. / 10., 7. / 240.,    0.};

    constexpr std::array<D, 11> weights4 = {
        0.,          7. / 240.,  -2. / 5., 169. / 60., -122. / 15., 91. / 8.,
        -122. / 15., 169. / 60., -2. / 5., 7. / 240.,  0.};

    constexpr std::array<D, 11> weights5 = {
        -13. / 288., 19. / 36., -87. / 32., 13. / 2.,   -323. / 48., 0.,
        323. / 48.,  -13. / 2., 87. / 32.,  -19. / 36., 13. / 288.};

    constexpr std::array<D, 11> weights6 = {
        13. / 240., -19. / 24., 87. / 16., -39. / 2.,  323. / 8., -1023. / 20.,
        323. / 8.,  -39. / 2.,  87. / 16., -19. / 24., 13. / 240.};

    std::array<D, 6> result{0., 0., 0., 0., 0., 0.};
    result[0] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights1.begin(), 0.) /
                epsilon;

    result[1] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights2.begin(), 0.) /
                std::pow(epsilon, 2);

    result[2] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights3.begin(), 0.) /
                std::pow(epsilon, 3);

    result[3] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights4.begin(), 0.) /
                std::pow(epsilon, 4);

    result[4] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights5.begin(), 0.) /
                std::pow(epsilon, 5);

    result[5] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights6.begin(), 0.) /
                std::pow(epsilon, 6);

    return result;
}

#endif // FINITE_DIFFERENCES_HPP
