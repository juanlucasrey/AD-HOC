#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include <cstdint>
#include "dco_codegen_helper.hpp"

#ifdef _MSC_VER
#pragma warning( push )
// do not warn on constant controlling expression
#pragma warning(disable:4127)
#endif

template <typename primal_t, typename derivative_t=primal_t>
void adjoint(
                  std::vector<primal_t> const& x,
                  std::vector<derivative_t>  & dx,
[[maybe_unused]]  std::vector<primal_t> const& p,
                  std::vector<primal_t>      & y,
                  std::vector<derivative_t>  & dy,
                  dco::codegen_buffer_t & ) {

  // Check vector sizes.
  if (x.size() < 4) {
    throw std::runtime_error("Input vector size too small.");
  }
  y.resize(1);
  if (dx.size() < 4) {
    throw std::runtime_error("Derivatives input vector size too small.");
  }
  if (dy.size() < 1) {
    throw std::runtime_error("Derivatives output vector size too small.");
  }


  // Forward code:
  const primal_t v0 = (x[2] * sqrt(x[3]));
  derivative_t dv0 = 0;
  const primal_t v1 = ((log((x[0] / x[1])) / v0) + (v0 * 5.00000000000000000e-01));
  derivative_t dv1 = 0;
  const primal_t v2 = (v1 - v0);
  derivative_t dv2 = 0;
  const primal_t v3 = (5.00000000000000000e-01 * erfc((v2 * -7.07106781186547462e-01)));
  derivative_t dv3 = 0;
  const primal_t v4 = (5.00000000000000000e-01 * erfc((v1 * -7.07106781186547462e-01)));
  derivative_t dv4 = 0;
  const primal_t v5 = ((x[0] * v4) - (x[1] * v3));
  derivative_t dv5 = 0;

  // Assigning outputs:
  y[0] = v5;

  // Assigning output adjoints:
  dv5  += dy[0];

  // Reverse code:

  // const primal_t v5 = ((x[0] * v4) - (x[1] * v3));
  dx[0] += v4 * dv5;
  dv4 += x[0] * dv5;
  dx[1] += ((-1.0)) * (v3) * dv5;
  dv3 += ((-1.0)) * (x[1]) * dv5;

  // const primal_t v4 = (5.00000000000000000e-01 * erfc((v1 * -7.07106781186547462e-01)));
  dv1 += (5.00000000000000000e-01) * ((-2.0) / sqrt((3.14159265358979323846264338327950288)) * exp(-(v1 * -7.07106781186547462e-01) * (v1 * -7.07106781186547462e-01))) * (-7.07106781186547462e-01) * dv4;

  // const primal_t v3 = (5.00000000000000000e-01 * erfc((v2 * -7.07106781186547462e-01)));
  dv2 += (5.00000000000000000e-01) * ((-2.0) / sqrt((3.14159265358979323846264338327950288)) * exp(-(v2 * -7.07106781186547462e-01) * (v2 * -7.07106781186547462e-01))) * (-7.07106781186547462e-01) * dv3;

  // const primal_t v2 = (v1 - v0);
  dv1 += dv2;
  dv0 += ((-1.0)) * dv2;

  // const primal_t v1 = ((log((x[0] / x[1])) / v0) + (v0 * 5.00000000000000000e-01));
  dx[0] += (1.0) / v0 * ((1.0) / (x[0] / x[1])) * ((1.0) / x[1]) * dv1;
  dx[1] += (1.0) / v0 * ((1.0) / (x[0] / x[1])) * (-(x[0] / x[1]) / x[1]) * dv1;
  dv0 += -(log((x[0] / x[1])) / v0) / v0 * dv1;
  dv0 += 5.00000000000000000e-01 * dv1;

  // const primal_t v0 = (x[2] * sqrt(x[3]));
  dx[2] += (sqrt(x[3])) * dv0;
  dx[3] += (x[2]) * ((1.0) / ((2.0) * sqrt(x[3]))) * dv0;

}

#ifdef _MSC_VC
#pragma warning( pop )
#endif
