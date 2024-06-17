
#include <base.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <combinatorics/combinations.hpp>
#include <combinatorics/factorial.hpp>
#include <combinatorics/integer_partition_index_sequence.hpp>
#include <combinatorics/partition_function.hpp>
#include <combinatorics/pow.hpp>
#include <utils/index_sequence.hpp>

#include <adhoc.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <complex>
#include <fstream>
#include <random>

namespace adhoc4 {

namespace {
std::array<double, 6> finite_differences(double x, double epsilon,
                                         std::function<double(double)> func) {
    std::vector<double> fvals{
        func(x - 5. * epsilon), func(x - 4. * epsilon), func(x - 3. * epsilon),
        func(x - 2. * epsilon), func(x - 1. * epsilon), func(x),
        func(x + 1. * epsilon), func(x + 2. * epsilon), func(x + 3. * epsilon),
        func(x + 4. * epsilon), func(x + 5. * epsilon)};

    std::vector<double> weights1 = {0.,        1. / 280., -4. / 105., 1. / 5.,
                                    -4. / 5.,  0.,        4. / 5.,    -1. / 5.,
                                    4. / 105., -1. / 280, 0.};

    std::vector<double> weights2 = {0.,        -1. / 560.,  8. / 315., -1. / 5.,
                                    8. / 5.,   -205. / 72., 8. / 5.,   -1. / 5.,
                                    8. / 315., -1. / 560.,  0.};

    std::vector<double> weights3 = {
        0.,         -7. / 240.,  3. / 10.,  -169. / 120., 61. / 30., 0.,
        -61. / 30., 169. / 120., -3. / 10., 7. / 240.,    0.};

    std::vector<double> weights4 = {
        0.,          7. / 240.,  -2. / 5., 169. / 60., -122. / 15., 91. / 8.,
        -122. / 15., 169. / 60., -2. / 5., 7. / 240.,  0.};

    std::vector<double> weights5 = {
        -13. / 288., 19. / 36., -87. / 32., 13. / 2.,   -323. / 48., 0.,
        323. / 48.,  -13. / 2., 87. / 32.,  -19. / 36., 13. / 288.};

    std::vector<double> weights6 = {
        13. / 240., -19. / 24., 87. / 16., -39. / 2.,  323. / 8., -1023. / 20.,
        323. / 8.,  -39. / 2.,  87. / 16., -19. / 24., 13. / 240.};

    std::array<double, 6> result;
    result.fill(0);
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
} // namespace

#define MACHEP 1.11022302462515654042E-16

double inline zeta(double x) {
    static double A[] = {
        12.0,
        -720.0,
        30240.0,
        -1209600.0,
        47900160.0,
        -1.8924375803183791606e9, /*1.307674368e12/691*/
        7.47242496e10,
        -2.950130727918164224e12,  /*1.067062284288e16/3617*/
        1.1646782814350067249e14,  /*5.109094217170944e18/43867*/
        -4.5979787224074726105e15, /*8.028576626982912e20/174611*/
        1.8152105401943546773e17,  /*1.5511210043330985984e23/854513*/
        -7.1661652561756670113e18  /*1.6938241367317436694528e27/236364091*/
    };

    int i;
    double a, b, k, s, t, w;

    // if (x == 1.0) {
    //     return std::numeric_limits<double>::infinity();
    // }

    if (x < 1.0) {
        return 0.;
    }

    // if (q <= 0.0) {
    //     if (q == floor(q)) {
    //         return std::numeric_limits<double>::infinity();
    //     }
    //     if (x != floor(x)) {
    //         return 0.;
    //     }
    // }

    /* Euler-Maclaurin summation formula */
    /*
    if( x < 25.0 )
    */
    {
        /* Permit negative q but continue sum until n+q > +9 .
         * This case should be handled by a reflection formula.
         * If q<0 and x is an integer, there is a relation to
         * the polygamma function.
         */
        s = std::pow(1.0, -x);
        a = 1.0;
        i = 0;
        b = 0.0;
        while ((i < 9) || (a <= 9.0)) {
            i += 1;
            a += 1.0;
            b = std::pow(a, -x);
            s += b;
            if (fabs(b / s) < MACHEP) {
                return s;
            }
        }

        w = a;
        s += b * w / (x - 1.0);
        s -= 0.5 * b;
        a = 1.0;
        k = 0.0;
        for (i = 0; i < 12; i++) {
            a *= x + k;
            b /= w;
            t = a * b / A[i];
            s = s + t;
            t = fabs(t / s);
            if (t < MACHEP) {
                return s;
            }
            k += 1.0;
            a *= x + k;
            b /= w;
            k += 1.0;
        }
        return s;
    }
}

double inline zeta2(double x) {

    double res = 0.;
    std::size_t i = 1;
    double contribution = std::pow(static_cast<double>(i), -x);
    i++;
    res += contribution;
    contribution = std::pow(static_cast<double>(i), -x);
    i++;
    while (std::abs(contribution / res) >
               std::numeric_limits<double>::epsilon() * 2. &&
           i < 10000) {
        res += contribution;
        contribution = std::pow(static_cast<double>(i), -x);
        i++;
    }

    // for (std::size_t i = 1; i < 200; i++) {
    //     res += std::pow(static_cast<double>(i), -x);
    // }
    return res;
}

const long double LOWER_THRESHOLD = 1.0e-6;
const long double UPPER_BOUND = 1.0e+4;
const int MAXNUM = 100;

std::complex<long double> zeta3(const std::complex<long double> &s) {
    std::complex<long double> a_arr[MAXNUM + 1];
    std::complex<long double> half(0.5, 0.0);
    std::complex<long double> one(1.0, 0.0);
    std::complex<long double> two(2.0, 0.0);
    std::complex<long double> rev(-1.0, 0.0);
    std::complex<long double> sum(0.0, 0.0);
    std::complex<long double> prev(1.0e+20, 0.0);

    a_arr[0] =
        half /
        (one -
         std::pow(two, (one - s))); // initialize with a_0 = 0.5 / (1 - 2^(1-s))
    sum += a_arr[0];

    for (int n = 1; n <= MAXNUM; n++) {
        std::complex<long double> nCplx(n, 0.0); // complex index

        for (int k = 0; k < n; k++) {
            std::complex<long double> kCplx(k, 0.0); // complex index

            a_arr[k] *= half * (nCplx / (nCplx - kCplx));
            sum += a_arr[k];
        }

        a_arr[n] =
            (rev * a_arr[n - 1] * std::pow((nCplx / (nCplx + one)), s) / nCplx);
        sum += a_arr[n];

        if (std::abs(prev - sum) <
            LOWER_THRESHOLD) // If the difference is less than or equal to the
                             // threshold value, it is considered to be
                             // convergent and the calculation is terminated.
            break;

        if (std::abs(sum) >
            UPPER_BOUND) // doesn't work for large values, so it gets terminated
                         // when it exceeds UPPER_BOUND
            break;

        prev = sum;
    }

    return sum;
}

long double zeta4(long double s) {
    long double a_arr[MAXNUM + 1];
    long double sum = 0.0;
    long double prev = 1.0e+20;

    // initialize with a_0 = 0.5 / (1 - 2^(1-s))
    a_arr[0] = 0.5 / (1.0 - std::pow(2.0, (1.0 - s)));
    sum += a_arr[0];

    // int n = 1;
    // for (int n = 1; n <= MAXNUM; n++) {
    for (int n = 1; n <= MAXNUM; n++) {
        long double nCplx(n); // index

        for (int k = 0; k < n; k++) {
            long double kCplx(k); // index

            a_arr[k] *= 0.5 * (nCplx / (nCplx - kCplx));
            sum += a_arr[k];
        }

        a_arr[n] = (-1.0 * a_arr[n - 1] * std::pow((nCplx / (nCplx + 1.0)), s) /
                    nCplx);
        sum += a_arr[n];

        // If the difference is less than or equal to the
        // threshold value, it is considered to be
        // convergent and the calculation is terminated.
        if (std::abs(prev - sum) < LOWER_THRESHOLD) {
            break;
        }

        // doesn't work for large values, so it gets terminated
        // when it exceeds UPPER_BOUND
        if (std::abs(sum) > UPPER_BOUND) {
            break;
        }

        prev = sum;
    }

    return sum;
}

long double zeta4d(long double s) {
    long double a_arr[MAXNUM + 1];
    long double a_arrd[MAXNUM + 1];
    long double sum = 0.0;
    long double prev = 1.0e+20;

    double const powval = std::pow(2.0, (1.0 - s));
    double const denom_sqrt = (1.0 - powval);
    a_arr[0] = 0.5 / denom_sqrt;
    a_arrd[0] =
        -0.5 * std::numbers::ln2_v<double> * powval / (denom_sqrt * denom_sqrt);
    sum += a_arrd[0];

    // int n = 1;
    for (int n = 1; n <= MAXNUM; n++) {
        long double nCplx(n); // index

        for (int k = 0; k < n; k++) {
            long double kCplx(k); // index

            a_arr[k] *= 0.5 * (nCplx / (nCplx - kCplx));
            a_arrd[k] *= 0.5 * (nCplx / (nCplx - kCplx));
            sum += a_arrd[k];
        }

        a_arr[n] = (-1.0 * a_arr[n - 1] * std::pow((nCplx / (nCplx + 1.0)), s) /
                    nCplx);

        a_arrd[n] = (-1.0 * a_arrd[n - 1] *
                     std::pow((nCplx / (nCplx + 1.0)), s) / nCplx);

        double const ln_temp = std::log(nCplx / (nCplx + 1.0));
        a_arrd[n] += (-1.0 * a_arr[n - 1] * ln_temp *
                      std::pow((nCplx / (nCplx + 1.0)), s) / nCplx);

        sum += a_arrd[n];

        // If the difference is less than or equal to the
        // threshold value, it is considered to be
        // convergent and the calculation is terminated.
        if (std::abs(prev - sum) < LOWER_THRESHOLD) {
            break;
        }

        // doesn't work for large values, so it gets terminated
        // when it exceeds UPPER_BOUND
        if (std::abs(sum) > UPPER_BOUND) {
            break;
        }

        prev = sum;
    }

    return sum;
}

template <std::size_t End, std::size_t CoeffIdx = 0, std::size_t Order,
          std::size_t... I>
void zeta4dd_aux3(std::index_sequence<I...> current, double &res,
                  std::array<double, Order> const &buff1,
                  std::array<double, Order> const &buff2) {
    constexpr auto sum_v = sum(current);
    constexpr auto coeff = BellCoeff(current);
    if constexpr (coeff == 1) {
        res += buff2[sum_v - 1] * inner_product_power(current, buff1);

    } else {
        res += buff2[sum_v - 1] * BellCoeff(current) *
               inner_product_power(current, buff1);
    }

    if constexpr ((CoeffIdx + 1) < End) {
        auto constexpr next = NextPartitionIS(current);
        zeta4dd_aux3<End, CoeffIdx + 1>(next, res, buff1, buff2);
    }
}

template <std::size_t Idx = 0, std::size_t Size>
void zeta4dd_aux2(double &res, std::array<double, Size> const &buff1,
                  std::array<double, Size> const &buff2) {
    res = 0;
    constexpr auto first = FirstPartitionIS<Idx + 1>();
    zeta4dd_aux3<partition_function(Idx + 1)>(first, res, buff1, buff2);
}

template <std::size_t Idx = 0, std::size_t Order, std::size_t Output>
void zeta4dd_aux(std::array<double, Output> &res,
                 std::array<double, Order> const &buff1,
                 std::array<double, Order> const &buff2) {
    zeta4dd_aux2<Idx>(res[Idx], buff1, buff2);
    if constexpr ((Idx + 1) < Order) {
        zeta4dd_aux<Idx + 1>(res, buff1, buff2);
    }
}

template <std::size_t End, std::size_t CoeffIdx = 0, std::size_t Order>
void mult_aux(double &res, double rightval,
              std::array<double, Order> const &right, double leftval,
              std::array<double, Order> const &left) {
    if constexpr (CoeffIdx == 0) {
        res += rightval * left[End - 1];
    } else if constexpr (CoeffIdx == End) {
        res += leftval * right[End - 1];
    } else {
        constexpr auto binomial = binomial_coefficient(End, CoeffIdx);
        res += binomial * left[CoeffIdx - 1] * right[End - 1 - CoeffIdx];
    }
    if constexpr (CoeffIdx < End) {
        mult_aux<End, CoeffIdx + 1>(res, rightval, right, leftval, left);
    }
}

template <std::size_t Idx = 0, std::size_t Order>
void mult(std::array<double, Order> &res, double rightval,
          std::array<double, Order> const &right, double leftval,
          std::array<double, Order> const &left) {
    res[Idx] = 0;
    mult_aux<Idx + 1>(res[Idx], rightval, right, leftval, left);
    if constexpr ((Idx + 1) < Order) {
        mult<Idx + 1>(res, rightval, right, leftval, left);
    }
}

template <std::size_t Order, std::size_t Output>
void zeta4dd(long double s, std::array<double, Output> &res) {
    long double a_arr[MAXNUM + 1];
    long double sum = 0.0;
    long double prev = 1.0e+20;

    std::array<double, Order> buffer;
    double const powval = std::pow(2.0, (1.0 - s));

    buffer[0] = powval * std::numbers::ln2_v<double>;
    for (std::size_t k = 1; k < Order; k++) {
        buffer[k] = buffer[k - 1] * -std::numbers::ln2_v<double>;
    }

    double const frac = 1.0 / (1.0 - powval);

    std::array<double, Order> buffer2;
    inv_t<double>::d<Order>(frac, frac, buffer2);

    zeta4dd_aux(res, buffer, buffer2);
    for (std::size_t k = 0; k < Order; k++) {
        res[k] *= 0.5;
    }

    a_arr[0] = 0.5 * frac;
    sum += a_arr[0];

    std::array<std::array<double, Order>, MAXNUM + 1> a_arrd;

    for (std::size_t k = 0; k < Order; k++) {
        a_arrd[0][k] = res[k];
    }

    // std::size_t n = 1;
    for (std::size_t n = 1; n <= MAXNUM; n++) {
        // for (int n = 1; n <= MAXNUM; n++) {
        long double nCplx(static_cast<long double>(n)); // index

        for (std::size_t k = 0; k < n; k++) {
            long double kCplx(static_cast<long double>(k)); // index
            double coeff = 0.5 * (nCplx / (nCplx - kCplx));

            a_arr[k] *= coeff;

            for (std::size_t d = 0; d < Order; d++) {
                a_arrd[k][d] *= coeff;
            }

            sum += a_arr[k];

            for (std::size_t d = 0; d < Order; d++) {
                res[d] += a_arrd[k][d];
            }
        }

        double const powval2 = std::pow((nCplx / (nCplx + 1.0)), s);
        double const coeff2 = -1.0 / nCplx;
        a_arr[n] = (-1.0 * a_arr[n - 1] * powval2 / nCplx);
        // a_arr[n] = a_arr[n - 1] * powval2;
        sum += a_arr[n];

        std::array<double, Order> buffer3;
        double const lnval = std::log(nCplx / (nCplx + 1.0));

        buffer3[0] = powval2 * lnval;
        for (std::size_t k = 1; k < Order; k++) {
            buffer3[k] = buffer3[k - 1] * lnval;
        }

        mult(a_arrd[n], powval2, buffer3, a_arr[n - 1], a_arrd[n - 1]);

        for (std::size_t d = 0; d < Order; d++) {
            a_arrd[n][d] *= coeff2;
            res[d] += a_arrd[n][d];
        }

        // If the difference is less than or equal to the
        // threshold value, it is considered to be
        // convergent and the calculation is terminated.
        // if (std::abs(prev - sum) < LOWER_THRESHOLD) {
        //     break;
        // }

        // doesn't work for large values, so it gets terminated
        // when it exceeds UPPER_BOUND
        // if (std::abs(sum) > UPPER_BOUND) {
        //     break;
        // }

        prev = sum;
    }
}

template <std::size_t Order, std::size_t Output>
static inline void zetad(double /* thisv */, double in,
                         std::array<double, Output> &res) {
    static_assert(Order <= Output);

    zeta4dd<Order>(in, res);
}

TEST(UnivariateFunctions, ZettaTimingOverOne) {
    std::uniform_real_distribution<double> d(1.1, 10.0);
    auto starttime = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    std::size_t iters = 100000;

    {
        std::mt19937 g(123);
        double acc = 0.;
        starttime = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < iters; i++) {
            double randomin = d(g);
            acc += std::riemann_zeta(randomin);
        }
        end = std::chrono::steady_clock::now();
        std::cout << "std "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
        acc /= static_cast<double>(iters);
        std::cout << acc << std::endl;
    }

    {
        std::mt19937 g(123);
        double acc = 0.;
        starttime = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < iters; i++) {
            double randomin = d(g);
            acc += zeta(randomin);
        }
        end = std::chrono::steady_clock::now();
        std::cout << "cephes "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
        acc /= static_cast<double>(iters);
        std::cout << acc << std::endl;
    }

    {
        std::mt19937 g(123);
        double acc = 0.;
        starttime = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < iters; i++) {
            double randomin = d(g);
            acc += zeta4(randomin);
        }
        end = std::chrono::steady_clock::now();
        std::cout << "stackoverflow "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
        acc /= static_cast<double>(iters);
        std::cout << acc << std::endl;
    }
}

TEST(UnivariateFunctions, ZettaTimingUnderOne) {
    std::uniform_real_distribution<double> d(-10, 0.9);
    auto starttime = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    std::size_t iters = 100000;

    {
        std::mt19937 g(123);
        double acc = 0.;
        starttime = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < iters; i++) {
            double randomin = d(g);
            acc += std::riemann_zeta(randomin);
        }
        end = std::chrono::steady_clock::now();
        std::cout << "std "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
        acc /= static_cast<double>(iters);
        std::cout << acc << std::endl;
    }

    // {
    //     std::mt19937 g(123);
    //     double acc = 0.;
    //     starttime = std::chrono::steady_clock::now();
    //     for (std::size_t i = 0; i < iters; i++) {
    //         double randomin = d(g);
    //         acc += zeta(randomin);
    //     }
    //     end = std::chrono::steady_clock::now();
    //     std::cout << "cephes "
    //               << std::chrono::duration_cast<std::chrono::microseconds>(
    //                      end - starttime)
    //                      .count()
    //               << " µs" << std::endl;
    //     acc /= static_cast<double>(iters);
    //     std::cout << acc << std::endl;
    // }

    {
        std::mt19937 g(123);
        double acc = 0.;
        starttime = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < iters; i++) {
            double randomin = d(g);
            acc += zeta4(randomin);
        }
        end = std::chrono::steady_clock::now();
        std::cout << "stackoverflow "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         end - starttime)
                         .count()
                  << " µs" << std::endl;
        acc /= static_cast<double>(iters);
        std::cout << acc << std::endl;
    }
}

inline auto zeta_borwein_positive(double x, std::size_t n) -> double {
    std::vector<double> d(n + 1, 0.);
    double prev = 0;
    for (std::size_t i = 0; i < (n + 1); i++) {
        // double numerator =
        //     std::pow(4, i) * factorial(static_cast<unsigned>(n + i - 1));
        // double denominateur = factorial(static_cast<unsigned>(n - i)) *
        //                       factorial(static_cast<unsigned>(2 * i));
        auto combinations = binomial_coefficient(n + i, 2 * i);
        // double temp1 = numerator / denominateur;
        double temp2 =
            static_cast<double>(combinations) / static_cast<double>(n + i);
        // double coeff = std::pow(4, i) *
        //                factorial(static_cast<unsigned>(n + i - 1)) /
        //                (factorial(static_cast<unsigned>(n - i)) *
        //                 factorial(static_cast<unsigned>(2 * i)));

        double coeff2 = std::pow(4, i) * temp2;
        d[i] = prev + coeff2;
        prev = d[i];
    }

    // for (std::size_t i = 0; i <= n; i++) {
    //     std::cout << d[i] << ", ";
    // }
    // std::cout << std::endl;

    double result = 0;
    for (std::size_t i = 0; i < n; i++) {
        double const_coeff = std::pow(-1., i) * (d[i] - d[n]);
        double const_coeff2 = const_coeff / d[n];
        double coeff = const_coeff2 * std::pow(i + 1, -x);
        result += coeff;
    }

    // result *= -1. / (d[n] * (1 - std::pow(2, 1 - x)));
    result *= -1. / (1 - std::pow(2, 1 - x));

    return result;
}

inline auto zeta_borwein(double x, std::size_t n) -> double {
    if (x < 0.) {
        double oneminuss = 1. - x;

        double sinpi = std::sin(0.5 * std::numbers::pi_v<double> * x);
        double pow = std::pow(2. * std::numbers::pi_v<double>, -oneminuss);
        double tgamm = std::tgamma(oneminuss);
        double zetatemp = zeta_borwein_positive(oneminuss, n);

        double result = sinpi * 2.0 * pow * tgamm * zetatemp;
        return result;
    }

    return zeta_borwein_positive(x, n);
}

TEST(UnivariateFunctions, ZetaPrecision) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    std::vector<double> in;
    std::vector<std::array<double, 7>> out;
    in.reserve(10000);
    out.reserve(10000);

    std::ifstream valuesfile;
    valuesfile.open("data/zeta.csv");
    ASSERT_TRUE(valuesfile.is_open());
    std::string line;
    std::string value;
    while (getline(valuesfile, line)) {
        std::istringstream tokenStream(line);
        std::getline(tokenStream, value, ',');
        in.push_back(std::stod(value));

        std::array<double, 7> results;

        std::getline(tokenStream, value, ',');
        results[0] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[1] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[2] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[3] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[4] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[5] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[6] = std::stod(value);

        out.push_back(results);
    }
    valuesfile.close();

    auto print_statistics = [&in, &out](auto &function) {
        double rel_average = 0.;
        double abs_rel_average = 0.;
        double rel_max = 0.;
        double abs_rel_max = 0.;
        std::size_t count = 0;
        for (std::size_t i = 0; i < in.size(); i++) {
            double const reference_out = out[i][0];
            double const reference_in = in[i];
            double const candidate_out = function(reference_in);
            double const error = candidate_out - reference_out;

            // double denom = reference_out;
            // if (denom == 0) {
            //     denom = candidate_out;
            // }

            if (reference_out != 0) {
                double const rel_error = error / reference_out;
                double const abs_rel_error = std::abs(rel_error);
                rel_average += rel_error;
                abs_rel_average += abs_rel_error;
                rel_max = std::max(rel_max, rel_error);
                abs_rel_max = std::max(abs_rel_max, abs_rel_error);
                count++;
            }
        }

        rel_average /= static_cast<double>(count);
        abs_rel_average /= static_cast<double>(count);

        std::cout << rel_average << std::endl;
        std::cout << abs_rel_average << std::endl;
        std::cout << rel_max << std::endl;
        std::cout << abs_rel_max << std::endl;
    };

    for (std::size_t i = 3; i < 25; i++) {
        std::cout << i << std::endl;
        std::function<double(double)> function = [i](double d) {
            return zeta_borwein(d, i);
        };
        print_statistics(function);
    }

    std::cout << "boost" << std::endl;
    std::function<double(double)> functionboost = [](double d) {
        return boost::math::zeta(d);
    };
    print_statistics(functionboost);

    std::cout << "std" << std::endl;
    std::function<double(double)> functionstd = [](double d) {
        return std::riemann_zeta(d);
    };
    print_statistics(functionstd);
}

template <std::size_t n>
inline auto zeta_borwein_positive_der(double x) -> double {
    std::vector<double> d(n + 1, 0.);
    double prev = 0;
    for (std::size_t i = 0; i < (n + 1); i++) {
        auto combinations = binomial_coefficient(n + i, 2 * i);
        double temp =
            static_cast<double>(combinations) / static_cast<double>(n + i);

        double coeff2 = std::pow(4, i) * temp;
        d[i] = prev + coeff2;
        prev = d[i];
    }

    double result = 0;
    for (std::size_t i = 0; i < n; i++) {
        double const_coeff = std::pow(-1., i) * (d[i] - d[n]);
        double const_coeff2 = const_coeff / d[n];
        double coeff = const_coeff2 * std::pow(i + 1, -x);
        result += coeff;
    }

    result *= -1. / (1 - std::pow(2, 1 - x));

    return result;
}

template <std::size_t Order, std::size_t n, std::size_t Output>
inline auto zeta_borwein_positive_der(double x, std::array<double, Output> &res)
    -> double {
    std::vector<double> d(n + 1, 0.);
    double prev = 0;
    for (std::size_t i = 0; i < (n + 1); i++) {
        auto combinations = binomial_coefficient(n + i, 2 * i);
        double temp =
            static_cast<double>(combinations) / static_cast<double>(n + i);

        double coeff2 = std::pow(4, i) * temp;
        d[i] = prev + coeff2;
        prev = d[i];
    }

    // double result = 1;
    std::array<double, Order> res1;
    res1.fill(0);
    double result = 0;
    // res[0] = 0;
    for (std::size_t i = 0; i < n; i++) {
        double const_coeff = std::pow(-1., i) * (d[i] - d[n]);
        double const_coeff2 = const_coeff / d[n];
        // double coeff = const_coeff2 / std::pow(i + 1, x);
        double coeff = const_coeff2 * std::pow(i + 1, -x);
        result += coeff;

        if (i != 0) {
            double temp = coeff;
            for (std::size_t j = 0; j < Order; j++) {
                temp *= -std::log(i + 1);
                res1[j] += temp;
            }
        }
    }

    std::array<double, Order> buffer;
    double const powval = std::pow(2.0, (1.0 - x));

    buffer[0] = powval * std::numbers::ln2_v<double>;
    for (std::size_t k = 1; k < Order; k++) {
        buffer[k] = buffer[k - 1] * -std::numbers::ln2_v<double>;
    }

    double const frac = 1.0 / (1.0 - powval);

    std::array<double, Order> buffer2;
    inv_t<double>::d<Order>(frac, frac, buffer2);

    std::array<double, Order> res2;
    zeta4dd_aux(res2, buffer, buffer2);

    for (std::size_t j = 0; j < Order; j++) {
        res2[j] *= -1.;
    }
    double result2 = -1. / (1 - std::pow(2, 1 - x));

    std::array<double, Order> resfinal;
    mult(resfinal, result2, res2, result, res1);
    result *= result2;

    for (std::size_t j = 0; j < Order; j++) {
        res[j] = resfinal[j];
    }

    return result;
}

TEST(UnivariateFunctions, ZetaBorweinDeriv) {

    std::array<double, 10> results1;
    results1.fill(0.);
    double val = 0.32;
    double res = zeta_borwein_positive_der<20>(val);

    zeta_borwein_positive_der<6, 20>(val, results1);
    // tgamma_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return zeta_borwein_positive_der<20>(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);

    EXPECT_LT(std::abs((results2[0] - results1[0]) / results2[0]), 1e-9);
    EXPECT_LT(std::abs((results2[1] - results1[1]) / results2[1]), 1e-9);
    EXPECT_LT(std::abs((results2[2] - results1[2]) / results2[2]), 1e-6);
    EXPECT_LT(std::abs((results2[3] - results1[3]) / results2[3]), 1e-6);
    EXPECT_LT(std::abs((results2[4] - results1[4]) / results2[4]), 1e-5);
    EXPECT_LT(std::abs((results2[5] - results1[5]) / results2[5]), 1e-5);
    EXPECT_EQ(results1[6], 0.);
    EXPECT_EQ(results1[7], 0.);
    EXPECT_EQ(results1[8], 0.);
    EXPECT_EQ(results1[9], 0.);
}

TEST(UnivariateFunctions, ZetaPrecisionDer) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    std::vector<double> in;
    std::vector<std::array<double, 7>> out;
    in.reserve(10000);
    out.reserve(10000);

    std::ifstream valuesfile;
    valuesfile.open("data/zeta.csv");
    ASSERT_TRUE(valuesfile.is_open());
    std::string line;
    std::string value;
    while (getline(valuesfile, line)) {
        std::istringstream tokenStream(line);
        std::getline(tokenStream, value, ',');
        in.push_back(std::stod(value));

        std::array<double, 7> results;

        std::getline(tokenStream, value, ',');
        results[0] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[1] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[2] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[3] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[4] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[5] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[6] = std::stod(value);

        out.push_back(results);
    }
    valuesfile.close();

    auto print_statistics = [&in, &out](auto &function, std::size_t n) {
        double rel_average = 0.;
        double abs_rel_average = 0.;
        double rel_max = 0.;
        double abs_rel_max = 0.;
        std::size_t count = 0;
        for (std::size_t i = 0; i < in.size(); i++) {
            double const reference_out = out[i][n];
            double const reference_in = in[i];
            double const candidate_out = function(reference_in, n);
            double const error = candidate_out - reference_out;

            // double denom = reference_out;
            // if (denom == 0) {
            //     denom = candidate_out;
            // }

            if (reference_out != 0) {
                double const rel_error = error / reference_out;
                double const abs_rel_error = std::abs(rel_error);
                rel_average += rel_error;
                abs_rel_average += abs_rel_error;
                rel_max = std::max(rel_max, rel_error);
                abs_rel_max = std::max(abs_rel_max, abs_rel_error);
                count++;
            }
        }

        rel_average /= static_cast<double>(count);
        abs_rel_average /= static_cast<double>(count);

        std::cout << rel_average << std::endl;
        std::cout << abs_rel_average << std::endl;
        std::cout << rel_max << std::endl;
        std::cout << abs_rel_max << std::endl;
    };

    std::function<double(double, std::size_t)> function = [](double d,
                                                             std::size_t i) {
        std::array<double, 10> results1;
        results1.fill(0.);
        zeta_borwein_positive_der<6, 20>(d, results1);

        return results1[i - 1];
    };

    for (std::size_t i = 1; i < 5; i++) {
        print_statistics(function, i);
    }
}

template <std::size_t n> inline auto zeta_borwein2(double x) -> double {
    if (x < 0.) {
        double oneminuss = 1. - x;

        double sinpi = std::sin(0.5 * std::numbers::pi_v<double> * x);
        double pow = std::pow(2. * std::numbers::pi_v<double>, -oneminuss);
        double tgamm = std::tgamma(oneminuss);
        double zetatemp = zeta_borwein_positive_der<n>(oneminuss);

        double result = sinpi * 2.0 * pow * tgamm * zetatemp;
        // double result = sinpi * pow * tgamm * zetatemp;
        return result;
    }

    return zeta_borwein_positive_der<n>(x);
}

template <std::size_t Order, std::size_t n, std::size_t Output>
inline auto zeta_borwein2_der(double x, std::array<double, Output> &res)
    -> double {
    if (x < 0.) {
        double oneminuss = 1. - x;

        constexpr double coeff = 0.5 * std::numbers::pi_v<double>;
        // double val = coeff * x;
        double xrem = std::remainder(x, 4);
        double val = coeff * xrem;
        double sinpi = std::sin(val);

        // double res = sin_t<double>::v(val);
        std::array<double, Order> sinders;
        sin_t<double>::d<Order>(sinpi, val, sinders);
        double multcoeff = coeff;
        for (std::size_t j = 0; j < Order; j++) {
            sinders[j] *= multcoeff;
            multcoeff *= coeff;
        }

        double pow = std::pow(2. * std::numbers::pi_v<double>, -oneminuss);

        std::array<double, Order> powders;
        const double ln2pi = std::log(2.0 * std::numbers::pi_v<double>);

        double multpow = ln2pi * pow;
        for (std::size_t j = 0; j < Order; j++) {
            powders[j] = multpow;
            multpow *= ln2pi;
        }

        std::array<double, Order> gammaders;
        double tgamm = std::tgamma(oneminuss);
        tgamma_t<double>::d<Order>(tgamm, oneminuss, gammaders);
        for (std::size_t j = 0; j < Order; j++) {
            gammaders[j] *= std::pow(-1., j + 1);
        }

        std::array<double, Order> zetaders;
        double zetatemp = zeta_borwein_positive_der<n>(oneminuss);
        zeta_borwein_positive_der<Order, n>(oneminuss, zetaders);
        for (std::size_t j = 0; j < Order; j++) {
            zetaders[j] *= std::pow(-1., j + 1);
        }

        std::array<double, Order> temp1;
        double val1 = tgamm * zetatemp;
        mult(temp1, zetatemp, zetaders, tgamm, gammaders);

        std::array<double, Order> temp2;
        double val2 = pow * val1;
        mult(temp2, val1, temp1, pow, powders);

        std::array<double, Order> temp3;
        double val3 = sinpi * val2;
        mult(temp3, val2, temp2, sinpi, sinders);

        for (std::size_t j = 0; j < Order; j++) {
            temp3[j] *= 2.0;
        }

        // double result = sinpi * 2.0 * pow * tgamm * zetatemp;
        double result = val3;
        for (std::size_t j = 0; j < Order; j++) {
            // res[j] = sinders[j];
            // res[j] = powders[j];
            res[j] = temp3[j];
        }
        return result;
    }

    return zeta_borwein_positive_der<Order, n>(x, res);
}

TEST(UnivariateFunctions, ZetaBorweinDerivNegative) {

    std::array<double, 10> results1;
    results1.fill(0.);
    double val = -0.32;
    double res = zeta_borwein2<20>(val);

    zeta_borwein2_der<6, 20>(val, results1);
    // tgamma_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return zeta_borwein2<20>(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);

    EXPECT_LT(std::abs((results2[0] - results1[0]) / results2[0]), 1e-9);
    EXPECT_LT(std::abs((results2[1] - results1[1]) / results2[1]), 1e-9);
    EXPECT_LT(std::abs((results2[2] - results1[2]) / results2[2]), 1e-6);
    EXPECT_LT(std::abs((results2[3] - results1[3]) / results2[3]), 1e-6);
    EXPECT_LT(std::abs((results2[4] - results1[4]) / results2[4]), 1e-5);
    EXPECT_LT(std::abs((results2[5] - results1[5]) / results2[5]), 1e-3);
    EXPECT_EQ(results1[6], 0.);
    EXPECT_EQ(results1[7], 0.);
    EXPECT_EQ(results1[8], 0.);
    EXPECT_EQ(results1[9], 0.);
}

TEST(UnivariateFunctions, ZetaPrecisionDerNegative) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    std::vector<double> in;
    std::vector<std::array<double, 7>> out;
    in.reserve(10000);
    out.reserve(10000);

    std::ifstream valuesfile;
    valuesfile.open("data/zeta.csv");
    ASSERT_TRUE(valuesfile.is_open());
    std::string line;
    std::string value;
    while (getline(valuesfile, line)) {
        std::istringstream tokenStream(line);
        std::getline(tokenStream, value, ',');
        in.push_back(std::stod(value));

        std::array<double, 7> results;

        std::getline(tokenStream, value, ',');
        results[0] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[1] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[2] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[3] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[4] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[5] = std::stod(value);
        std::getline(tokenStream, value, ',');
        results[6] = std::stod(value);

        out.push_back(results);
    }
    valuesfile.close();

    auto print_statistics = [&in, &out](auto &function, std::size_t n) {
        double rel_average = 0.;
        double abs_rel_average = 0.;
        double rel_max = 0.;
        double abs_rel_max = 0.;
        std::size_t count = 0;
        for (std::size_t i = 0; i < in.size(); i++) {
            double const reference_out = out[i][n];
            double const reference_in = in[i];
            double const candidate_out = function(reference_in, n);

            // double temp = zeta_borwein2<20>(reference_in);

            double const error = candidate_out - reference_out;

            // double denom = reference_out;
            // if (denom == 0) {
            //     denom = candidate_out;
            // }

            if (reference_out != 0) {
                double const rel_error = error / reference_out;
                double const abs_rel_error = std::abs(rel_error);
                rel_average += rel_error;
                abs_rel_average += abs_rel_error;
                rel_max = std::max(rel_max, rel_error);
                abs_rel_max = std::max(abs_rel_max, abs_rel_error);
                count++;
            }
        }

        rel_average /= static_cast<double>(count);
        abs_rel_average /= static_cast<double>(count);

        std::cout << rel_average << std::endl;
        std::cout << abs_rel_average << std::endl;
        std::cout << rel_max << std::endl;
        std::cout << abs_rel_max << std::endl;
    };

    std::function<double(double, std::size_t)> function = [](double d,
                                                             std::size_t i) {
        std::array<double, 10> results1;
        results1.fill(0.);
        zeta_borwein2_der<6, 20>(d, results1);

        return results1[i - 1];
    };

    for (std::size_t i = 1; i < 5; i++) {
        print_statistics(function, i);
    }
}

} // namespace adhoc4
