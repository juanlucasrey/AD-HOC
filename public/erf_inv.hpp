/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/adhoc2).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ADHOC_ERF_INV_HPP
#define ADHOC_ERF_INV_HPP

#include "constants_constexpr.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace adhoc {

enum class algorithm {
    bsm,
    acklam,
    hybrid,
    voultier_narrow,
    voultier_wide,
    boost
};

namespace detail {

template <bool Precise, auto Method, bool CheckInput> struct quantile {};

template <bool Precise, auto Method> struct quantile_sn_half_domain {};

// J. D. Beasley and S. G. Springer, The percentage points of the Normal
// Distribution, Applied Statistics .26(1977), 118–121.
inline auto bsm_core(double x) -> double {
    constexpr double a0 = 2.50662823884;
    constexpr double a1 = -18.61500062529;
    constexpr double a2 = 41.39119773534;
    constexpr double a3 = -25.44106049637;
    constexpr double b0 = -8.47351093090;
    constexpr double b1 = 23.08336743743;
    constexpr double b2 = -21.06224101826;
    constexpr double b3 = 3.13082909833;
    x = 0.5 - x;
    const double t = x * x;
    return x * (((a3 * t + a2) * t + a1) * t + a0) /
           ((((b3 * t + b2) * t + b1) * t + b0) * t + 1.0);
}

inline auto bsm_tail(double x) -> double {
    constexpr double a0 = 0.3374754822726147;
    constexpr double a1 = 0.9761690190917186;
    constexpr double a2 = 0.1607979714918209;
    constexpr double a3 = 0.0276438810333863;
    constexpr double a4 = 0.0038405729373609;
    constexpr double a5 = 0.0003951896511919;
    constexpr double a6 = 0.0000321767881768;
    constexpr double a7 = 0.0000002888167364;
    constexpr double a8 = 0.0000003960315187;
    x = std::log(-std::log(x));
    return (((((((a8 * x + a7) * x + a6) * x + a5) * x + a4) * x + a3) * x +
             a2) *
                x +
            a1) *
               x +
           a0;
}

// https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
inline auto acklam_core(double x) -> double {
    constexpr double a0 = 2.506628277459239;
    constexpr double a1 = -30.66479806614716;
    constexpr double a2 = 138.3577518672690;
    constexpr double a3 = -275.9285104469687;
    constexpr double a4 = 220.9460984245205;
    constexpr double a5 = -39.69683028665376;
    constexpr double b0 = -13.28068155288572;
    constexpr double b1 = 66.80131188771972;
    constexpr double b2 = -155.6989798598866;
    constexpr double b3 = 161.5858368580409;
    constexpr double b4 = -54.47609879822406;

    x = 0.5 - x;
    const double t = x * x;
    return x * (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t + a0) /
           (((((b4 * t + b3) * t + b2) * t + b1) * t + b0) * t + 1.0);
}

inline auto acklam_tail(double x) -> double {
    constexpr double a0 = -377.42823772627105;
    constexpr double a1 = -561.95698649904261;
    constexpr double a2 = 327.53143278065073;
    constexpr double a3 = 308.39462027700063;
    constexpr double a4 = 41.414137442352711;
    constexpr double a5 = 1.0000254722080659;
    constexpr double b0 = 128.45717255699017;
    constexpr double b1 = 482.28072133209952;
    constexpr double b2 = 314.09501777996508;
    constexpr double b3 = 41.423215642907323;

    x = std::sqrt(-2.0 * std::log(x));
    return (((((a5 * x + a4) * x + a3) * x + a2) * x + a1) * x + a0) /
           ((((x + b3) * x + b2) * x + b1) * x + b0);
}

// Paul M. Voutier, A New Approximation to the Normal Distribution Quantile
// Function, https://arxiv.org/pdf/1002.0567.pdf, 2010-02-03
inline auto voutier_core_narrow(double x) -> double {
    constexpr double a0 = 0.195740115269792;
    constexpr double a1 = -0.652871358365296;
    constexpr double a2 = 1.246899760652504;
    constexpr double b0 = 0.155331081623168;
    constexpr double b1 = -0.839293158122257;
    x = 0.5 - x;
    const double t = x * x;
    return x * (a2 + (a1 * t + a0) / ((t + b1) * t + b0));
}

inline auto voutier_core_wide(double x) -> double {
    constexpr double a0 = 0.151015505647689;
    constexpr double a1 = -0.5303572634357367;
    constexpr double a2 = 1.365020122861334;
    constexpr double b0 = 0.132089632343748;
    constexpr double b1 = -0.7607324991323768;
    x = 0.5 - x;
    const double t = x * x;
    return x * (a2 + (a1 * t + a0) / ((t + b1) * t + b0));
}

inline auto voutier_tail(double x) -> double {
    constexpr double a0 = -16.682320830719986527;
    constexpr double a1 = -4.120411523939115059;
    constexpr double a2 = -0.029814187308200211;
    constexpr double a3 = 1.000182518730158122;
    constexpr double b0 = 7.173787663925508066;
    constexpr double b1 = 8.759693508958633869;
    x = std::sqrt(-2.0 * std::log(x));
    return a3 * x + a2 + (a1 * x + a0) / ((x + b1) * x + b0);
}

inline auto boost_core(double x) -> double {
    constexpr double a0 = -0.81190817642945057;
    constexpr double a1 = -13.354748714112132;
    constexpr double a2 = 53.428042578062737;
    constexpr double a3 = -20.254723493920086;
    constexpr double a4 = -58.348072057523751;
    constexpr double a5 = 35.087978088960163;
    constexpr double a6 = 13.128355134221572;
    constexpr double a7 = -8.5976748167713133;
    constexpr double a8 = 0.12605094006504738;
    constexpr double b0 = 1128.3925808289032;
    constexpr double b1 = -1094.5464942300521;
    constexpr double b2 = -1766.7756985800709;
    constexpr double b3 = 1762.7924746230976;
    constexpr double b4 = 747.36694965761842;
    constexpr double b5 = -803.74164942783432;
    constexpr double b6 = -59.511016498857209;
    constexpr double b7 = 89.739221245048256;
    constexpr double b8 = -2.6335978648960228;
    x = 1.0 - 2.0 * x;
    return x * (x + 10.0) *
           (a8 +
            (((((((a7 * x + a6) * x + a5) * x + a4) * x + a3) * x + a2) * x +
              a1) *
                 x +
             a0) /
                (((((((((x + b8) * x + b7) * x + b6) * x + b5) * x + b4) * x +
                    b3) *
                       x +
                   b2) *
                      x +
                  b1) *
                     x +
                 b0));
}

inline auto boost_middle(double x) -> double {
    constexpr double a0 = 0.13333711696293529;
    constexpr double a1 = 9.8865163242534972;
    constexpr double a2 = 114.49467335699001;
    constexpr double a3 = 152.77499813600113;
    constexpr double a4 = -579.07865566125281;
    constexpr double a5 = 10.960477450244705;
    constexpr double a6 = 535.27764392638426;
    constexpr double a7 = -261.16446549436466;
    constexpr double a8 = 17.231572625270488;
    constexpr double b0 = 0.10780240715673521;
    constexpr double b1 = 8.4688484971192519;
    constexpr double b2 = 106.78266162252052;
    constexpr double b3 = 203.04273505253553;
    constexpr double b4 = -447.36263414682918;
    constexpr double b5 = -187.20979177931818;
    constexpr double b6 = 472.11174615018859;
    constexpr double b7 = -151.20949797583876;
    x *= 2.0;
    return std::sqrt(-std::log(x)) *
           (((((((((a8 * x + a7) * x + a6) * x + a5) * x + a4) * x + a3) * x +
               a2) *
                  x +
              a1) *
                 x +
             a0) /
            ((((((((x + b7) * x + b6) * x + b5) * x + b4) * x + b3) * x + b2) *
                  x +
              b1) *
                 x +
             b0));
}

inline auto boost_tail_1(double x) -> double {
    constexpr double a0 = 0.00018458909362031794;
    constexpr double a1 = 1.7261027424073558;
    constexpr double a2 = 4.7554115183418926;
    constexpr double a3 = 10.2449261556033;
    constexpr double a4 = 13.466420170635713;
    constexpr double a5 = 12.395358851586725;
    constexpr double a6 = 8.2736359631587515;
    constexpr double a7 = 1.4169150014184384;
    constexpr double a8 = -0.00012877767817177246;
    constexpr double a9 = 4.6272593862766304e-06;
    constexpr double a10 = -8.7102848588613201e-08;
    constexpr double b0 = 1.3794866041755256;
    constexpr double b1 = 3.7780225165695414;
    constexpr double b2 = 8.9366675004102216;
    constexpr double b3 = 12.409906778685142;
    constexpr double b4 = 14.262344541658194;
    constexpr double b5 = 10.398902654783781;
    constexpr double b6 = 5.8930624273244803;
    return x *
           (((((((((((a10 * x + a9) * x + a8) * x + a7) * x + a6) * x + a5) *
                     x +
                 a4) *
                    x +
                a3) *
                   x +
               a2) *
                  x +
              a1) *
                 x +
             a0) /
            (((((((x + b6) * x + b5) * x + b4) * x + b3) * x + b2) * x + b1) *
                 x +
             b0));
}

inline auto boost_tail_2(double x) -> double {
    constexpr double a0 = -0.17111535636706701;
    constexpr double a1 = 25.296324910372075;
    constexpr double a2 = 55.745341300112059;
    constexpr double a3 = 87.475838697719723;
    constexpr double a4 = 87.632073009946481;
    constexpr double a5 = 23.321790380202494;
    constexpr double a6 = 1.4145607038332064;
    constexpr double a7 = -5.4433574088557471e-06;
    constexpr double a8 = 4.9257580493180608e-08;
    constexpr double b0 = 19.625301967454732;
    constexpr double b1 = 44.81943246551964;
    constexpr double b2 = 82.415498497767231;
    constexpr double b3 = 83.266409024188391;
    constexpr double b4 = 64.115832746021297;
    constexpr double b5 = 16.506368822563267;
    return x *
           (((((((((a8 * x + a7) * x + a6) * x + a5) * x + a4) * x + a3) * x +
               a2) *
                  x +
              a1) *
                 x +
             a0) /
            ((((((x + b5) * x + b4) * x + b3) * x + b2) * x + b1) * x + b0));
}

inline auto boost_tail_3(double x) -> double {
    constexpr double a0 = -5294.7319526873998;
    constexpr double a1 = -7118.2883266390218;
    constexpr double a2 = -1814.7538057088468;
    constexpr double a3 = -39.12188527380232;
    constexpr double a4 = 19.891240779361034;
    constexpr double a5 = 1.4187058222621047;
    constexpr double a6 = 0.023185357020071816;
    constexpr double a7 = -1.6468682414767039e-07;
    constexpr double a8 = 4.9633061970674137e-10;
    constexpr double a9 = 1.3910604440247263;
    constexpr double b0 = 4260.1227576901019;
    constexpr double b1 = 7112.8220969169251;
    constexpr double b2 = 10805.590039994282;
    constexpr double b3 = 5779.1014677624871;
    constexpr double b4 = 1028.9591543532936;
    constexpr double b5 = 61.552572414050204;
    return x *
           (a9 +
            ((((((((a8 * x + a7) * x + a6) * x + a5) * x + a4) * x + a3) * x +
               a2) *
                  x +
              a1) *
                 x +
             a0) /
                ((((((x + b5) * x + b4) * x + b3) * x + b2) * x + b1) * x +
                 b0));
}

inline auto boost_tail_4(double x) -> double {
    constexpr double a0 = -69889.333605522406;
    constexpr double a1 = 375426.85685930907;
    constexpr double a2 = 590951.23596874031;
    constexpr double a3 = 151980.50529542295;
    constexpr double a4 = 11033.491028195533;
    constexpr double a5 = 247.97166208835344;
    constexpr double a6 = 1.4142165427779521;
    constexpr double a7 = -3.2393724616396713e-09;
    constexpr double b0 = 232231.02997271836;
    constexpr double b1 = 406539.91047666734;
    constexpr double b2 = 432611.80731446017;
    constexpr double b3 = 107903.54942394678;
    constexpr double b4 = 7805.1783691859673;
    constexpr double b5 = 175.34391657125514;
    return x *
           ((((((((a7 * x + a6) * x + a5) * x + a4) * x + a3) * x + a2) * x +
              a1) *
                 x +
             a0) /
            ((((((x + b5) * x + b4) * x + b3) * x + b2) * x + b1) * x + b0));
}

template <> struct quantile_sn_half_domain<false, algorithm::bsm> {
    static inline auto impl(double x) -> double {
        if (x > 0.08) {
            return bsm_core(x);
        }

        return bsm_tail(x);
    }
};

template <> struct quantile_sn_half_domain<false, algorithm::acklam> {
    static inline auto impl(double x) -> double {
        if (x >= 0.02425) {
            return acklam_core(x);
        }

        return acklam_tail(x);
    }
};

template <> struct quantile_sn_half_domain<false, algorithm::hybrid> {
    static inline auto impl(double x) -> double {
        if (x > 0.08) {
            return bsm_core(x);
        }

        return quantile_sn_half_domain<false, algorithm::acklam>::impl(x);
    }
};

template <> struct quantile_sn_half_domain<false, algorithm::voultier_narrow> {
    static inline auto impl(double x) -> double {
        if (x > 0.0465) {
            return voutier_core_narrow(x);
        }

        return voutier_tail(x);
    }
};

template <> struct quantile_sn_half_domain<false, algorithm::voultier_wide> {
    static inline auto impl(double x) -> double {
        if (x > 0.025) {
            return voutier_core_wide(x);
        }

        return voutier_tail(x);
    }
};

template <> struct quantile_sn_half_domain<false, algorithm::boost> {
    static inline auto impl(double x) -> double {
        if (x >= 0.25) {
            return boost_core(x);
        }

        if (x >= 0.125) {
            return boost_middle(x);
        }

        x = std::sqrt(-std::log(2.0 * x));

        if (x < 3) {
            return boost_tail_1(x);
        }

        if (x < 6) {
            return boost_tail_2(x);
        }

        if (x < 18) {
            return boost_tail_3(x);
        }

        return boost_tail_4(x);
    }
};

template <class D = double> inline auto cdf(D const &x) {
    using constexpression::sqrt;
    using std::erfc; // to allow Argument-dependent lookup
    constexpr double one_over_root_two = 1.0 / sqrt(2.0);
    return 0.5 * erfc(x * -one_over_root_two);
};

template <auto Method> struct quantile_sn_half_domain<true, Method> {
    static inline auto impl(double x) -> double {
        // below this small value, halley corection overflows and underflows,
        // giving NaN
        if (x >= 4.3458473798968777e-311) { // 2^(-1031)
            const double r = -quantile_sn_half_domain<false, Method>::impl(x);
            using constexpression::pi;
            using constexpression::sqrt;
            constexpr double root_two_pi = sqrt(2.0 * pi<double>());
            x = (cdf(r) - x) * root_two_pi * std::exp(r * r * 0.5);
            return -r + x / (1.0 + r * x * 0.5); // Halley's method
        }

        // boost method is the best we got at this range
        // (0, 4.3458473798968777e-311)
        return quantile_sn_half_domain<false, algorithm::boost>::impl(x);
    }
};

template <bool Precise, auto Method> struct quantile<Precise, Method, false> {
    static inline auto impl(double x) -> double {
        const double sign = std::copysign(1.0, x - 0.5);
        x = sign > 0. ? 1.0 - x : x;
        return quantile_sn_half_domain<Precise, Method>::impl(x) * sign;
    }
};

template <bool Precise, auto Method> struct quantile<Precise, Method, true> {
    static inline auto impl(double x) -> double {
        if (x > 1.0 || x < 0.0) {
            throw std::invalid_argument("quantile input must be in [0,1]");
        }

        if (x == 1.0) {
            return std::numeric_limits<double>::infinity();
        }

        if (x == 0.0) {
            return -std::numeric_limits<double>::infinity();
        }

        return quantile<Precise, Method, false>::impl(x);
    }
};

template <auto Method, class ForwardIt> struct quantile_vec {
    static inline void impl(ForwardIt first, ForwardIt last) {
        // default implementation
        while (first != last) {
            detail::quantile<false, Method, false>::impl(*first++);
        }
    }
};

template <class ForwardIt> struct quantile_vec<algorithm::bsm, ForwardIt> {
    static inline void impl(ForwardIt first, ForwardIt last) {
        constexpr double a0 = 2.50662823884;
        constexpr double a1 = -18.61500062529;
        constexpr double a2 = 41.39119773534;
        constexpr double a3 = -25.44106049637;
        constexpr double b0 = -8.47351093090;
        constexpr double b1 = 23.08336743743;
        constexpr double b2 = -21.06224101826;
        constexpr double b3 = 3.13082909833;

        // we vectorise core calculations
        auto it = first;
        while (it != last) {
            double &x = *it++;
            x -= 0.5;
            const double sign = std::copysign(1.0, x);
            x = 0.5 - std::abs(x);

            const double separation = std::copysign(1.0, x - 0.08);
            const double is_core = (separation + 1.0) * 0.5;
            const double is_tail = (1.0 - separation) * 0.5;

            x = is_core * 0.5 -
                separation * x; // to maintain x value intact if tail value

            const double t = x * x;
            x = x * is_tail +
                x * is_core * (((a3 * t + a2) * t + a1) * t + a0) /
                    ((((b3 * t + b2) * t + b1) * t + b0) * t + 1.0);

            // core results cannot go beyond 1.98 so we add 2.0 to tail input
            // values to make sure we recognise them in the next loop
            const double add = 2.0 * is_tail;
            x += add;
            x *= sign;
        }

        // tail calculations cannot be vectorised because they use sqrt and log
        while (first != last) {
            double &x = *first++;
            if (std::abs(x) > 1.99) {
                const double sign = std::copysign(1.0, x);
                x = std::abs(x) - 2.0;
                x = bsm_tail(x) * sign;
            }
        }
    }
};

template <class ForwardIt> struct quantile_vec<algorithm::acklam, ForwardIt> {
    static inline void impl(ForwardIt first, ForwardIt last) {
        constexpr double a0 = 2.506628277459239;
        constexpr double a1 = -30.66479806614716;
        constexpr double a2 = 138.3577518672690;
        constexpr double a3 = -275.9285104469687;
        constexpr double a4 = 220.9460984245205;
        constexpr double a5 = -39.69683028665376;
        constexpr double b0 = -13.28068155288572;
        constexpr double b1 = 66.80131188771972;
        constexpr double b2 = -155.6989798598866;
        constexpr double b3 = 161.5858368580409;
        constexpr double b4 = -54.47609879822406;

        // we vectorise core calculations
        auto it = first;
        while (it != last) {
            double &x = *it++;
            x -= 0.5;
            const double sign = std::copysign(1.0, x);
            x = 0.5 - std::abs(x);

            const double separation = std::copysign(1.0, x - 0.02425);
            const double is_core = (separation + 1.0) * 0.5;
            const double is_tail = (1.0 - separation) * 0.5;

            x = is_core * 0.5 -
                separation * x; // to maintain x value intact if tail value

            const double t = x * x;
            x = x * is_tail +
                x * is_core *
                    (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t + a0) /
                    (((((b4 * t + b3) * t + b2) * t + b1) * t + b0) * t + 1.0);

            // core results cannot go beyond 1.98 so we add 2.0 to tail input
            // values to make sure we recognise them in the next loop
            const double add = 2.0 * is_tail;
            x += add;
            x *= sign;
        }

        // tail calculations cannot be vectorised because they use sqrt and log
        while (first != last) {
            double &x = *first++;
            if (std::abs(x) > 1.99) {
                const double sign = std::copysign(1.0, x);
                x = std::abs(x) - 2.0;
                x = acklam_tail(x) * sign;
            }
        }
    }
};

} // namespace detail

template <bool Precise = false, auto Method = algorithm::bsm,
          bool CheckInput = false>
inline auto erfc_inv(double x) -> double {
    using constexpression::sqrt;
    constexpr double minus_one_over_root_two = -1.0 / sqrt(2.0);
    return minus_one_over_root_two *
           detail::quantile<Precise, Method, CheckInput>::impl(x * 0.5);
}

template <bool Precise = false, auto Method = algorithm::bsm,
          bool CheckInput = false>
inline auto erf_inv(double x) -> double {
    using constexpression::sqrt;
    constexpr double one_over_root_two = 1.0 / sqrt(2.0);
    return one_over_root_two *
           detail::quantile<Precise, Method, CheckInput>::impl((x + 1.0) * 0.5);
}

template <auto Method, class ForwardIt>
inline void erfc_inv(ForwardIt first, ForwardIt last) {
    auto it = first;
    while (it != last) {
        *it++ *= 0.5;
    }
    detail::quantile_vec<Method, ForwardIt>::impl(first, last);

    using constexpression::sqrt;
    constexpr double minus_one_over_root_two = -1.0 / sqrt(2.0);

    it = first;
    while (it != last) {
        *it++ *= minus_one_over_root_two;
    }
}

template <auto Method, class ForwardIt>
inline void erf_inv(ForwardIt first, ForwardIt last) {
    auto it = first;
    while (it != last) {
        *it++ = (*it + 1.0) * 0.5;
    }
    detail::quantile_vec<Method, ForwardIt>::impl(first, last);

    using constexpression::sqrt;
    constexpr double one_over_root_two = 1.0 / sqrt(2.0);

    it = first;
    while (it != last) {
        *it++ *= one_over_root_two;
    }
}

} // namespace adhoc

#endif // ADHOC_ERF_INV_HPP
