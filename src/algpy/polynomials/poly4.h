#ifndef ALGPY_POLY4_H
#define ALGPY_POLY4_H

#define ALGPY_POLY4_N 4

#include <float.h>

#include "poly3.h"

/* Ascending real coefficients; roots_to_coeffs returns the monic polynomial.
 * Its input roots must be real or occur in conjugate pairs, in any order.
 */
static inline void
poly4_roots_to_coeffs(const complex result[ALGPY_POLY4_N], double p[ALGPY_POLY4_N+1])
{
    double ar = result[0].data[0], ai = result[0].data[1];
    double br = result[1].data[0], bi = result[1].data[1];
    double cr = result[2].data[0], ci = result[2].data[1];
    double dr = result[3].data[0], di = result[3].data[1];
    double abr = ar*br - ai*bi, abi = ar*bi + ai*br;
    double cdr = cr*dr - ci*di, cdi = cr*di + ci*dr;
    double sr = ar+br, si = ai+bi, tr = cr+dr, ti = ci+di;
    p[0] = abr*cdr - abi*cdi;
    p[1] = -sr*cdr + si*cdi - tr*abr + ti*abi;
    p[2] = abr + cdr + sr*tr - si*ti;
    p[3] = -sr-tr;
    p[4] = 1;
}

/* p[4] must be nonzero. Ferrari factorization uses one resolvent cubic and
 * two quadratics. This is a small fixed-cost helper for algebra eigenvalues. */
static inline void
poly4_coeffs_to_roots(const double p[ALGPY_POLY4_N+1], complex result[ALGPY_POLY4_N])
{
    double inv = 1/p[4];
    double a = p[3]*inv, b = p[2]*inv, c = p[1]*inv, d = p[0]*inv;
    double shift = a/4;
    double pp = b - 6*shift*shift;
    double qq = (8*shift*shift - 2*b)*shift + c;
    double rr = ((-3*shift*shift + b)*shift - c)*shift + d;
    double qscale = fabs(8*shift*shift*shift) + fabs(2*b*shift) + fabs(c);
    /* In a centred polynomial, qq may itself consist only of accumulated
     * multiplication roundoff. Include the natural cubic coefficient scale;
     * otherwise a nearly repeated biquadratic can give s=0 and qq/s below. */
    qscale += fabs(pp)*sqrt(fabs(pp)) + sqrt(fabs(rr))*sqrt(sqrt(fabs(rr)));
    int i;
    if (d == 0) {
        poly3_coeffs_to_roots(p+1, result);
        result[3] = complex_from_re(0);
        return;
    }
    if (fabs(qq) <= 4*DBL_EPSILON*qscale) {
        /* Biquadratic in y = x + a/4, allowing roundoff in qq. Degenerate
         * Clifford algebras can produce this case with repeated eigenvalues. */
        double quadratic[3] = {rr, pp, 1};
        complex squares[2], u, v;
        poly2_coeffs_to_roots(quadratic, squares);
        u = complex_sqrt(squares[0]);
        v = complex_sqrt(squares[1]);
        result[0] = complex_make(u.data[0]-shift, u.data[1]);
        result[1] = complex_make(-u.data[0]-shift, -u.data[1]);
        result[2] = complex_make(v.data[0]-shift, v.data[1]);
        result[3] = complex_make(-v.data[0]-shift, -v.data[1]);
    }
    else {
        /* m = s^2: m^3 + 2*pp*m^2 + (pp^2 - 4*rr)*m - qq^2 = 0. */
        double cubic[4] = {-qq*qq, pp*pp - 4*rr, 2*pp, 1};
        complex roots[3];
        double m, s, left[3], right[3];
        poly3_coeffs_to_roots(cubic, roots);
        m = roots[0].data[0];
        for (i = 1; i < 3; i++) {
            if (roots[i].data[1] == 0 && roots[i].data[0] > m) {
                m = roots[i].data[0];
            }
        }
        s = sqrt(m);
        left[0] = (pp+m-qq/s)/2;
        left[1] = s;
        left[2] = 1;
        right[0] = (pp+m+qq/s)/2;
        right[1] = -s;
        right[2] = 1;
        poly2_coeffs_to_roots(left, result);
        poly2_coeffs_to_roots(right, result+2);
        for (i = 0; i < 4; i++) {
            result[i].data[0] -= shift;
        }
    }
}

#endif /* ALGPY_POLY4_H */
