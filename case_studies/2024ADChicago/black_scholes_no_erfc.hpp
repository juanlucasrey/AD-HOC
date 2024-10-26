#ifndef CASE_STUDIES_2024ADCHICAGO_BLACK_SCHOLES_NO_ERFC_HPP
#define CASE_STUDIES_2024ADCHICAGO_BLACK_SCHOLES_NO_ERFC_HPP

#include <cmath>
#include <numbers>

template <class D> inline auto cdf_n(D const &x) {
    using std::erf;
    constexpr double minus_one_over_root_two =
        -1.0 / std::numbers::sqrt2_v<double>;
    return 0.5 * (1.0 - erf(x * minus_one_over_root_two));
}

template <class I1, class I2, class I3, class I4>
inline auto call_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * 0.5;
    auto d2 = d1 - totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

template <class I1, class I2, class I3, class I4>
inline auto put_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * 0.5;
    auto d2 = d1 - totalvol;
    return K * cdf_n(-d2) - S * cdf_n(-d1);
}

#endif // CASE_STUDIES_2024ADCHICAGO_BLACK_SCHOLES_NO_ERFC_HPP
