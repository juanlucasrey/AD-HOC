#ifndef TEST_CALL_PRICE4_HPP
#define TEST_CALL_PRICE4_HPP

#include <constants_type.hpp>

#include <cmath>
#include <numbers>

namespace adhoc4 {

template <class D> inline auto cdf_n(D const &x) {
    using constants::CD;
    using std::erfc;
    constexpr double minus_one_over_root_two =
        -1.0 / std::numbers::sqrt2_v<double>;
    return CD<0.5>() * erfc(x * CD<minus_one_over_root_two>());
}

template <class I1, class I2, class I3, class I4>
auto call_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<0.5>();
    auto d2 = d1 - totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

template <class I1, class I2, class I3, class I4>
auto put_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<0.5>();
    auto d2 = d1 - totalvol;
    return K * cdf_n(-d2) - S * cdf_n(-d1);
}

} // namespace adhoc4

#endif // TEST_CALL_PRICE4_HPP
