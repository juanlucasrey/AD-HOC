#ifndef CASE_STUDIES_2024ADCHICAGO_BLACK_SCHOLES_HPP
#define CASE_STUDIES_2024ADCHICAGO_BLACK_SCHOLES_HPP

#ifdef WITHOUT_REF
#define DECL_INPUT(D) D
#else
#define DECL_INPUT(D) D const &
#endif

#ifdef WITHOUT_AUTO
#define DECL_COND_AUTO(D) D
#else
#define DECL_COND_AUTO(D) auto
#endif

#include "../../include/constants_type.hpp"

#include <cmath>
#include <numbers>

template <class D> inline auto cdf_n(DECL_INPUT(D) x) {
    using adhoc::constants::CD;
    using std::erfc;
    constexpr double minus_one_over_root_two =
        -1.0 / std::numbers::sqrt2_v<double>;
    return CD<0.5>() * erfc(x * CD<minus_one_over_root_two>());
}

template <class I1, class I2, class I3, class I4>
inline auto call_price(DECL_INPUT(I1) S, DECL_INPUT(I2) K, DECL_INPUT(I3) v,
                       DECL_INPUT(I4) T) {
    using adhoc::constants::CD;
    using std::log;
    using std::sqrt;
    DECL_COND_AUTO(I1) totalvol = v * sqrt(T);
    DECL_COND_AUTO(I1) d1 = log(S / K) / totalvol + totalvol * CD<0.5>();
    auto d2 = d1 - totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

template <class I1, class I2, class I3, class I4>
inline auto put_price(DECL_INPUT(I1) S, DECL_INPUT(I2) K, DECL_INPUT(I3) v,
                      DECL_INPUT(I4) T) {
    using adhoc::constants::CD;
    using std::log;
    using std::sqrt;
    DECL_COND_AUTO(I1) totalvol = v * sqrt(T);
    DECL_COND_AUTO(I1) d1 = log(S / K) / totalvol + totalvol * CD<0.5>();
    auto d2 = d1 - totalvol;
    return K * cdf_n(-d2) - S * cdf_n(-d1);
}

#endif // CASE_STUDIES_2024ADCHICAGO_BLACK_SCHOLES_HPP
