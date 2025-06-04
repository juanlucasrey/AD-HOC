
/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef BRNG_DISTRIBUTION_NORMAL_DISTRIBUTION
#define BRNG_DISTRIBUTION_NORMAL_DISTRIBUTION

#include "centered_canonical_distribution.hpp"

#include <cmath>
#include <limits>
#include <type_traits>

namespace adhoc {

namespace detail {

template <class RealType = double>
inline auto boost_core(RealType x) -> RealType {
    constexpr RealType a0 = -0.81190817642945057;
    constexpr RealType a1 = -13.354748714112132;
    constexpr RealType a2 = 53.428042578062737;
    constexpr RealType a3 = -20.254723493920086;
    constexpr RealType a4 = -58.348072057523751;
    constexpr RealType a5 = 35.087978088960163;
    constexpr RealType a6 = 13.128355134221572;
    constexpr RealType a7 = -8.5976748167713133;
    constexpr RealType a8 = 0.12605094006504738;
    constexpr RealType b0 = 1128.3925808289032;
    constexpr RealType b1 = -1094.5464942300521;
    constexpr RealType b2 = -1766.7756985800709;
    constexpr RealType b3 = 1762.7924746230976;
    constexpr RealType b4 = 747.36694965761842;
    constexpr RealType b5 = -803.74164942783432;
    constexpr RealType b6 = -59.511016498857209;
    constexpr RealType b7 = 89.739221245048256;
    constexpr RealType b8 = -2.6335978648960228;
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

template <class RealType = double>
inline auto boost_middle(RealType x) -> RealType {
    constexpr RealType a0 = 0.13333711696293529;
    constexpr RealType a1 = 9.8865163242534972;
    constexpr RealType a2 = 114.49467335699001;
    constexpr RealType a3 = 152.77499813600113;
    constexpr RealType a4 = -579.07865566125281;
    constexpr RealType a5 = 10.960477450244705;
    constexpr RealType a6 = 535.27764392638426;
    constexpr RealType a7 = -261.16446549436466;
    constexpr RealType a8 = 17.231572625270488;
    constexpr RealType b0 = 0.10780240715673521;
    constexpr RealType b1 = 8.4688484971192519;
    constexpr RealType b2 = 106.78266162252052;
    constexpr RealType b3 = 203.04273505253553;
    constexpr RealType b4 = -447.36263414682918;
    constexpr RealType b5 = -187.20979177931818;
    constexpr RealType b6 = 472.11174615018859;
    constexpr RealType b7 = -151.20949797583876;
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

template <class RealType = double>
inline auto boost_tail_1(RealType x) -> RealType {
    constexpr RealType a0 = 0.00018458909362031794;
    constexpr RealType a1 = 1.7261027424073558;
    constexpr RealType a2 = 4.7554115183418926;
    constexpr RealType a3 = 10.2449261556033;
    constexpr RealType a4 = 13.466420170635713;
    constexpr RealType a5 = 12.395358851586725;
    constexpr RealType a6 = 8.2736359631587515;
    constexpr RealType a7 = 1.4169150014184384;
    constexpr RealType a8 = -0.00012877767817177246;
    constexpr RealType a9 = 4.6272593862766304e-06;
    constexpr RealType a10 = -8.7102848588613201e-08;
    constexpr RealType b0 = 1.3794866041755256;
    constexpr RealType b1 = 3.7780225165695414;
    constexpr RealType b2 = 8.9366675004102216;
    constexpr RealType b3 = 12.409906778685142;
    constexpr RealType b4 = 14.262344541658194;
    constexpr RealType b5 = 10.398902654783781;
    constexpr RealType b6 = 5.8930624273244803;
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

template <class RealType = double>
inline auto boost_tail_2(RealType x) -> RealType {
    constexpr RealType a0 = -0.17111535636706701;
    constexpr RealType a1 = 25.296324910372075;
    constexpr RealType a2 = 55.745341300112059;
    constexpr RealType a3 = 87.475838697719723;
    constexpr RealType a4 = 87.632073009946481;
    constexpr RealType a5 = 23.321790380202494;
    constexpr RealType a6 = 1.4145607038332064;
    constexpr RealType a7 = -5.4433574088557471e-06;
    constexpr RealType a8 = 4.9257580493180608e-08;
    constexpr RealType b0 = 19.625301967454732;
    constexpr RealType b1 = 44.81943246551964;
    constexpr RealType b2 = 82.415498497767231;
    constexpr RealType b3 = 83.266409024188391;
    constexpr RealType b4 = 64.115832746021297;
    constexpr RealType b5 = 16.506368822563267;
    return x *
           (((((((((a8 * x + a7) * x + a6) * x + a5) * x + a4) * x + a3) * x +
               a2) *
                  x +
              a1) *
                 x +
             a0) /
            ((((((x + b5) * x + b4) * x + b3) * x + b2) * x + b1) * x + b0));
}

template <class RealType = double>
inline auto boost_tail_3(RealType x) -> RealType {
    constexpr RealType a0 = -5294.7319526873998;
    constexpr RealType a1 = -7118.2883266390218;
    constexpr RealType a2 = -1814.7538057088468;
    constexpr RealType a3 = -39.12188527380232;
    constexpr RealType a4 = 19.891240779361034;
    constexpr RealType a5 = 1.4187058222621047;
    constexpr RealType a6 = 0.023185357020071816;
    constexpr RealType a7 = -1.6468682414767039e-07;
    constexpr RealType a8 = 4.9633061970674137e-10;
    constexpr RealType a9 = 1.3910604440247263;
    constexpr RealType b0 = 4260.1227576901019;
    constexpr RealType b1 = 7112.8220969169251;
    constexpr RealType b2 = 10805.590039994282;
    constexpr RealType b3 = 5779.1014677624871;
    constexpr RealType b4 = 1028.9591543532936;
    constexpr RealType b5 = 61.552572414050204;
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

template <class RealType = double>
inline auto boost_tail_4(RealType x) -> RealType {
    constexpr RealType a0 = -69889.333605522406;
    constexpr RealType a1 = 375426.85685930907;
    constexpr RealType a2 = 590951.23596874031;
    constexpr RealType a3 = 151980.50529542295;
    constexpr RealType a4 = 11033.491028195533;
    constexpr RealType a5 = 247.97166208835344;
    constexpr RealType a6 = 1.4142165427779521;
    constexpr RealType a7 = -3.2393724616396713e-09;
    constexpr RealType b0 = 232231.02997271836;
    constexpr RealType b1 = 406539.91047666734;
    constexpr RealType b2 = 432611.80731446017;
    constexpr RealType b3 = 107903.54942394678;
    constexpr RealType b4 = 7805.1783691859673;
    constexpr RealType b5 = 175.34391657125514;
    return x *
           ((((((((a7 * x + a6) * x + a5) * x + a4) * x + a3) * x + a2) * x +
              a1) *
                 x +
             a0) /
            ((((((x + b5) * x + b4) * x + b3) * x + b2) * x + b1) * x + b0));
}

template <class RealType = double> inline auto impl(RealType x) -> RealType {
    if (x >= 0.25) {
        return boost_core<RealType>(x);
    }

    if (x >= 0.125) {
        return boost_middle<RealType>(x);
    }

    x = std::sqrt(-std::log(2.0 * x));

    if (x < 3) {
        return boost_tail_1<RealType>(x);
    }

    if (x < 6) {
        return boost_tail_2<RealType>(x);
    }

    if (x < 18) {
        return boost_tail_3<RealType>(x);
    }

    return boost_tail_4<RealType>(x);
}

} // namespace detail

template <class RealType = double> class standard_normal_distribution {

    static_assert(std::is_floating_point_v<RealType>);

  public:
    using result_type = RealType;

    standard_normal_distribution() = default;

    void reset() {}

    template <class RNG> auto operator()(RNG &g) -> result_type {
        adhoc::centered_canonical_distribution<RealType, false> u;
        auto const x = u(g);

        const RealType sign = std::copysign(1.0, -x);
        return detail::impl<RealType>(std::abs(x)) * sign;
    }

    static constexpr auto mean() -> result_type { return 0.; }
    static constexpr auto stddev() -> result_type { return 1.; }

    static constexpr auto min() -> result_type {
        return -std::numeric_limits<RealType>::max();
    }
    static constexpr auto max() -> result_type {
        return std::numeric_limits<RealType>::max();
    }

    constexpr auto
    operator==(const standard_normal_distribution & /* rhs */) const -> bool {
        return true;
    }

    constexpr auto
    operator!=(const standard_normal_distribution & /* rhs */) const -> bool {
        return false;
    }
};

template <class RealType = double> class normal_distribution {

    static_assert(std::is_floating_point_v<RealType>);

  public:
    using result_type = RealType;

    normal_distribution() : normal_distribution(0.0) {}

    explicit normal_distribution(RealType mean, RealType stddev = 1.0)
        : _mean(mean), _stddev(stddev) {}

    void reset() {}

    template <class RNG> auto operator()(RNG &g) -> result_type {
        standard_normal_distribution sn;
        return (sn(g) * _stddev) + _mean;
    }

    auto mean() const -> result_type { return _mean; }
    auto stddev() const -> result_type { return _stddev; }

    static constexpr auto min() -> result_type {
        return -std::numeric_limits<RealType>::max();
    }
    static constexpr auto max() -> result_type {
        return std::numeric_limits<RealType>::max();
    }

    auto operator==(const normal_distribution &rhs) const -> bool {
        return (this->_mean == rhs._mean) && (this->_stddev == rhs._stddev);
    }

    auto operator!=(const normal_distribution &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    result_type _mean;
    result_type _stddev;
};

template <class RNG, class RealType = double> class standard_normal {
  public:
    using result_type = RealType;

    explicit standard_normal() : rng_() {}

    template <class... Args>
    explicit standard_normal(Args const &...args) : rng_(args...) {}

    template <class... Args>
    explicit standard_normal(Args &...args) : rng_(args...) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        auto x = this->rng_.template operator()<FwdDirection>();
        const RealType sign = std::copysign(1.0, -x);
        return detail::impl<RealType>(std::abs(x)) * sign;
    };

    void discard(unsigned long long z) { rng_.discard(z); }

    auto operator==(standard_normal const &rhs) const noexcept -> bool {
        return this->rng_ == rhs.rng_;
    }

    auto operator!=(const standard_normal &rhs) const -> bool {
        return !(this->operator==(rhs));
    }

  private:
    centered_canonical_open<RNG, RealType> rng_{RNG()};
};

} // namespace adhoc

#endif // BRNG_DISTRIBUTION_NORMAL_DISTRIBUTION
