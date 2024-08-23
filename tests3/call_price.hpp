#ifndef CALL_PRICE_HPP
#define CALL_PRICE_HPP

#include <constants_constexpr.hpp>
#include <constants_type.hpp>

#include <cmath>

namespace adhoc3 {

template <class D> inline auto cdf_n(D const &x) {
    using constants::CD;
    using constants::encode;
    using std::erfc;
    constexpr double minus_one_over_root_two =
        -1.0 / constexpression::sqrt(2.0);
    return CD<encode(0.5)>() * erfc(x * CD<encode(minus_one_over_root_two)>());
    // return CD<0.5>() * erfc(x * CD<minus_one_over_root_two>());
}

template <class I1, class I2, class I3, class I4>
auto call_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using constants::encode;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<encode(0.5)>();
    auto d2 = d1 - totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

template <class I1, class I2, class I3, class I4>
auto digital_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using constants::encode;
    using std::exp;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d2 = log(S / K) / totalvol - totalvol * CD<encode(0.5)>();
    // auto d2 = d1 - totalvol;
    return cdf_n(d2);
}

template <class I1, class I2, class I3, class I4>
auto put_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using constants::encode;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<encode(0.5)>();
    auto d2 = d1 - totalvol;
    return K * cdf_n(-d2) - S * cdf_n(-d1);
}

} // namespace adhoc3

#endif // CALL_PRICE_HPP
