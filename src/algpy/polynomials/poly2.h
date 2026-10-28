#ifndef ALGPY_POLY2_H
#define ALGPY_POLY2_H

#define ALGPY_POLY2_N 2

#include "../fields/complex.h"

/* Ascending real coefficients; roots_to_coeffs returns the monic polynomial.
 * Its input roots must be real or occur in conjugate pairs, in any order.
 */

static inline void
poly2_roots_to_coeffs(const complex result[ALGPY_POLY2_N], double p[ALGPY_POLY2_N+1])
{
    double ar = result[0].data[0], ai = result[0].data[1];
    double br = result[1].data[0], bi = result[1].data[1];
    p[0] = ar*br - ai*bi;
    p[1] = -ar-br;
    p[2] = 1;
}

/* p[2] must be nonzero. Use the root product to avoid cancellation in
 * the smaller real root of the quadratic formula. */
static inline void
poly2_coeffs_to_roots(const double p[ALGPY_POLY2_N+1], complex result[ALGPY_POLY2_N])
{
    double b = p[1]/p[2], c = p[0]/p[2];
    double d = b*b - 4*c;
    if (d >= 0) {
        double q = -(b + copysign(sqrt(d), b))/2;
        result[0] = complex_make(q, 0);
        result[1] = complex_make(q == 0 ? 0 : c/q, 0);
    }
    else {
        double im = sqrt(-d)/2;
        result[0] = complex_make(-b/2, im);
        result[1] = complex_make(-b/2, -im);
    }
}

#endif /* ALGPY_POLY2_H */
