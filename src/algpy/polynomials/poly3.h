#ifndef ALGPY_POLY3_H
#define ALGPY_POLY3_H

#define ALGPY_POLY3_N 3

#include "poly2.h"

/* Ascending real coefficients; roots_to_coeffs returns the monic polynomial.
 * Its input roots must be real or occur in conjugate pairs, in any order.
 */

static inline void
poly3_roots_to_coeffs(const complex result[ALGPY_POLY3_N], double p[ALGPY_POLY3_N+1])
{
    double ar = result[0].data[0], ai = result[0].data[1];
    double br = result[1].data[0], bi = result[1].data[1];
    double cr = result[2].data[0], ci = result[2].data[1];
    double abr = ar*br - ai*bi, abi = ar*bi + ai*br;
    p[0] = -abr*cr + abi*ci;
    p[1] = abr + (ar+br)*cr - (ai+bi)*ci;
    p[2] = -ar-br-cr;
    p[3] = 1;
}

/* p[3] must be nonzero. Cardano for one real root, the trigonometric
 * formula for three real roots. All arithmetic stays in double precision. */
static inline void
poly3_coeffs_to_roots(const double p[ALGPY_POLY3_N+1], complex result[ALGPY_POLY3_N])
{
    double a = p[2]/p[3], b = p[1]/p[3], c = p[0]/p[3];
    double shift = a/3;
    double u = (b - a*shift)/3;
    double v = ((2*shift*shift - b)*shift + c)/2;
    double d = v*v + u*u*u;
    if (c == 0) {
        double quadratic[3] = {b, a, 1};
        result[0] = complex_from_re(0);
        poly2_coeffs_to_roots(quadratic, result+1);
    }
    else if (d > 0) {
        double r = cbrt(-v - copysign(sqrt(d), v));
        double s = r == 0 ? 0 : -u/r;
        double re = -(r+s)/2 - shift;
        double im = sqrt(3.0)*(r-s)/2;
        result[0] = complex_make(r+s-shift, 0);
        result[1] = complex_make(re, im);
        result[2] = complex_make(re, -im);
    }
    else if (d == 0) {
        double r = cbrt(-v);
        result[0] = complex_make(2*r-shift, 0);
        result[1] = result[2] = complex_make(-r-shift, 0);
    }
    else {
        double r = sqrt(-u);
        double theta = acos(fmax(-1, fmin(1, -v/(r*r*r))))/3;
        double x = r*cos(theta), y = sqrt(3.0)*r*sin(theta);
        result[0] = complex_make(2*x-shift, 0);
        result[1] = complex_make(-x-y-shift, 0);
        result[2] = complex_make(-x+y-shift, 0);
    }
}

#endif /* ALGPY_POLY3_H */
