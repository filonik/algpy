#ifndef ALGPY_POLY1_H
#define ALGPY_POLY1_H

#define ALGPY_POLY1_N 1

#include "../fields/complex.h"

/* Coefficients are in ascending powers. Roots must describe a real
 * polynomial (real roots and conjugate pairs); the result is monic.
 */

static inline void
poly1_roots_to_coeffs(const complex result[ALGPY_POLY1_N], double p[ALGPY_POLY1_N+1])
{
    p[0] = -result[0].data[0];
    p[1] = 1;
}

/* p[1] must be nonzero. */
static inline void
poly1_coeffs_to_roots(const double p[ALGPY_POLY1_N+1], complex result[ALGPY_POLY1_N])
{
    result[0] = complex_make(-p[0]/p[1], 0);
}

#endif /* ALGPY_POLY1_H */
