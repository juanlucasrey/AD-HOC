#ifndef BLACK_SCHOLES_C_INCLUDED
#define BLACK_SCHOLES_C_INCLUDED

#include <math.h>
#include <stdio.h>

double cdf_n(double x) { return 0.5 * erfc(x * -0.70710678118654746); }

double call_price(double S, double K, double v, double T) {
    double totalvol = v * sqrt(T);
    double d1 = log(S / K) / totalvol + totalvol * 0.5;
    double d2 = d1 - totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

#endif
