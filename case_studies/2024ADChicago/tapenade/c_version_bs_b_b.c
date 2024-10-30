/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.16 (feature_julia) - 28 Oct 2024 17:01
*/
// #include <adStack.h>
/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.16 (feature_julia) - 28 Oct 2024 17:01
*/
#include "adStack.h"
#include <math.h>
#include <stdio.h>

/*
  Differentiation of cdf_n_b in reverse (adjoint) mode:
   gradient     of useful results: *xb
   with respect to varying inputs: x cdf_nb
   Plus diff mem management of: xb:in


  Differentiation of cdf_n in reverse (adjoint) mode:
   gradient     of useful results: cdf_n
   with respect to varying inputs: x
*/
void cdf_n_b_b(double x, double *xb0, double *xb, double *xbb, double cdf_nb,
        double *cdf_nbb) {
    double temp;
    temp = -(0.70710678118654746*0.70710678118654746*(x*x));
    *xb0 = -(2*x*pow(0.70710678118654746, 3)*exp(temp)*cdf_nb*1.1283791671*0.5
        *(*xbb));
    *cdf_nbb = exp(temp)*1.1283791671*0.5*0.70710678118654746*(*xbb);
}

/*
  Differentiation of cdf_n in reverse (adjoint) mode:
   gradient     of useful results: cdf_n
   with respect to varying inputs: x
*/
void cdf_n_b_nodiff(double x, double *xb, double cdf_nb) {
    *xb = 0.70710678118654746*1.1283791671*exp(-(-(0.70710678118654746*x)*-(
        0.70710678118654746*x)))*0.5*cdf_nb;
}

/*
  Differentiation of cdf_n_nodiff in reverse (adjoint) mode:
   gradient     of useful results: x cdf_n_nodiff
   with respect to varying inputs: x
*/
void cdf_n_nodiff_b(double x, double *xb, double cdf_n_nodiffb) {
    *xb = *xb + 0.70710678118654746*1.1283791671*exp(-(-(0.70710678118654746*x
        )*-(0.70710678118654746*x)))*0.5*cdf_n_nodiffb;
}

double cdf_n_nodiff_nodiff(double x) {
    return 0.5*erfc(x*-0.70710678118654746);
}

/*
  Differentiation of call_price_b in reverse (adjoint) mode:
   gradient     of useful results: *Tb v *Kb *vb call_priceb *Sb
                K S T
   with respect to varying inputs: *Tb v *Kb *vb call_priceb *Sb
                K S T
   RW status of diff variables: Tb:(loc) *Tb:in-out v:incr Kb:(loc)
                *Kb:in-out vb:(loc) *vb:in-out call_priceb:incr
                Sb:(loc) *Sb:in-out K:incr S:incr T:incr
   Plus diff mem management of: Tb:in Kb:in vb:in Sb:in


  Differentiation of call_price in reverse (adjoint) mode:
   gradient     of useful results: v call_price K S T
   with respect to varying inputs: v K S T
   RW status of diff variables: v:incr call_price:in-killed K:incr
                S:incr T:incr
*/
void call_price_b_b(double S, double *Sb0, double *Sb, double *Sbb, double K,
        double *Kb0, double *Kb, double *Kbb, double v, double *vb0, double *
        vb, double *vbb, double T, double *Tb0, double *Tb, double *Tbb,
        double call_priceb, double *call_pricebb) {
    double totalvol = v*sqrt(T);
    double totalvolb0 = 0.0;
    double totalvolb = 0.0;
    double totalvolbb = 0.0;
    double d1 = log(S/K)/totalvol + totalvol*0.5;
    double d1b0 = 0.0;
    double d1b = 0.0;
    double d1bb = 0.0;
    double d2 = d1 - totalvol;
    double d2b0 = 0.0;
    double d2b = 0.0;
    double d2bb = 0.0;
    double result1;
    double result1b0;
    double result1b;
    double result1bb;
    double result2;
    double result2b0;
    double result2b;
    double result2bb;
    double temp;
    double tempb1;
    double tempb;
    double tempbb;
    double temp0;
    double tempb0;
    result1 = cdf_n_nodiff_nodiff(d1);
    result2 = cdf_n_nodiff_nodiff(d2);
    result1b = S*call_priceb;
    result2b = -(K*call_priceb);
    d2b = 0.0;
    cdf_n_b_nodiff(d2, &d2b, result2b);
    d1b = 0.0;
    cdf_n_b_nodiff(d1, &d1b, result1b);
    d1b = d1b + d2b;
    temp = S/K;
    tempb = d1b/(K*temp*totalvol);
    totalvolb = (0.5-log(temp)/(totalvol*totalvol))*d1b - d2b;
    temp = sqrt(T);
    if (!(T==0.0)) {
        tempb0 = *Tbb/(2.0*temp);
        *vb0 = *vb0 + totalvolb*tempb0;
        totalvolbb = v*tempb0;
        tempb1 = -(v*totalvolb*tempb0/temp);
    } else {
        tempb1 = 0.0;
        totalvolbb = 0.0;
    }
    tempb1 = tempb1 + totalvolb*(*vbb);
    totalvolbb = totalvolbb + temp*(*vbb);
    if (!(T==0.0))
        *Tb0 = *Tb0 + tempb1/(2.0*sqrt(T));
    temp = S/K;
    temp0 = log(temp)/(totalvol*totalvol);
    tempb0 = -(d1b*totalvolbb/(totalvol*totalvol));
    d1bb = (0.5-temp0)*totalvolbb;
    tempb1 = tempb0/temp;
    totalvolb0 = -(2*totalvol*temp0*tempb0);
    result2b0 = -(call_priceb*(*Kbb));
    tempbb = *Sbb - temp*(*Kbb);
    result1b0 = call_priceb*(*Sbb);
    temp0 = K*temp*totalvol;
    d1bb = d1bb + tempbb/temp0;
    d2bb = d1bb - totalvolbb;
    tempb0 = -(d1b*tempbb/(temp0*temp0));
    tempb1 = tempb1 + K*totalvol*tempb0 - tempb*(*Kbb);
    *Kb0 = *Kb0 + temp*totalvol*tempb0;
    d1b0 = 0.0;
    result1bb = 0.0;
    cdf_n_b_b(d1, &d1b0, &d1b, &d1bb, result1b, &result1bb);
    d2b0 = 0.0;
    result2bb = 0.0;
    cdf_n_b_b(d2, &d2b0, &d2b, &d2bb, result2b, &result2bb);
    *call_pricebb = *call_pricebb + result1*(*Sbb) - result2*(*Kbb) + S*
        result1bb - K*result2bb;
    cdf_n_nodiff_b(d2, &d2b0, result2b0);
    cdf_n_nodiff_b(d1, &d1b0, result1b0);
    d1b0 = d1b0 + d2b0;
    temp0 = S/K;
    totalvolb0 = totalvolb0 + K*temp*tempb0 + (0.5-log(temp0)/(totalvol*
        totalvol))*d1b0 - d2b0;
    tempb0 = d1b0/(K*temp0*totalvol);
    *Sb0 = *Sb0 + tempb1/K + call_priceb*result1bb + tempb0;
    *Kb0 = *Kb0 - S*tempb1/(K*K) - call_priceb*result2bb - temp0*tempb0;
    temp0 = sqrt(T);
    *vb0 = *vb0 + temp0*totalvolb0;
    if (!(T==0.0))
        *Tb0 = *Tb0 + v*totalvolb0/(2.0*temp0);
}
