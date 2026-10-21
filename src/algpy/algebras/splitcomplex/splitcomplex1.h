/*
 * Arithmetic of split complex numbers (x + y*j), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_SPLITCOMPLEX1_H
#define ALGPY_SPLITCOMPLEX1_H

#define ALGPY_SPLITCOMPLEX1_N 2

#include <math.h>

#include "../../vectorspaces/vec2.h"

typedef vec2 splitcomplex1;

static inline splitcomplex1
splitcomplex1_make(double x, double y)
{
    splitcomplex1 r = {{x, y}};
    return r;
}

static inline splitcomplex1
splitcomplex1_re(splitcomplex1 a)
{
    return splitcomplex1_make(a.data[0], 0);
}

static inline splitcomplex1
splitcomplex1_im(splitcomplex1 a)
{
    return splitcomplex1_make(0, a.data[1]);
}

static inline splitcomplex1
splitcomplex1_from_re(double a)
{
    return splitcomplex1_make(a, 0);
}

static inline double
splitcomplex1_to_re(splitcomplex1 a)
{
    return a.data[0];
}

static inline splitcomplex1
splitcomplex1_one(void)
{
    return splitcomplex1_from_re(1);
}

static inline splitcomplex1
splitcomplex1_mul(splitcomplex1 x, splitcomplex1 y)
{
    return splitcomplex1_make(
        x.data[0]*y.data[0] + x.data[1]*y.data[1],
        x.data[0]*y.data[1] + x.data[1]*y.data[0]);
}

static inline void
splitcomplex1_basis(splitcomplex1 result[ALGPY_SPLITCOMPLEX1_N])
{
    result[0] = splitcomplex1_make(1, 0);
    result[1] = splitcomplex1_make(0, 1);
}

static inline void
splitcomplex1_coeffs(
    splitcomplex1 result[ALGPY_SPLITCOMPLEX1_N * ALGPY_SPLITCOMPLEX1_N])
{
    splitcomplex1 basis[ALGPY_SPLITCOMPLEX1_N];
    int i, j;
    splitcomplex1_basis(basis);
    for (i = 0; i < ALGPY_SPLITCOMPLEX1_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX1_N; j++) {
            result[i * ALGPY_SPLITCOMPLEX1_N + j] =
                splitcomplex1_mul(basis[i], basis[j]);
        }
    }
}

static inline splitcomplex1
splitcomplex1_conjugate(splitcomplex1 a)
{
    return splitcomplex1_make(a.data[0], -a.data[1]);
}

/* Norm a*conj(a) = x^2 - y^2 */
static inline double
splitcomplex1_norm2(splitcomplex1 a)
{
    return a.data[0]*a.data[0] - a.data[1]*a.data[1];
}

static inline splitcomplex1
splitcomplex1_inv(splitcomplex1 a)
{
    double n = splitcomplex1_norm2(a);
    return splitcomplex1_make(a.data[0] / n, -a.data[1] / n);
}

static inline splitcomplex1
splitcomplex1_div(splitcomplex1 a, splitcomplex1 b)
{
    return splitcomplex1_mul(a, splitcomplex1_inv(b));
}

static inline splitcomplex1
splitcomplex1_exp(splitcomplex1 a)
{
    double e = exp(a.data[0]);
    return splitcomplex1_make(e * cosh(a.data[1]),
                              e * sinh(a.data[1]));
}

static inline splitcomplex1
splitcomplex1_log(splitcomplex1 a)
{
    double p = a.data[0] + a.data[1];
    double m = a.data[0] - a.data[1];
    double lp = log(p);
    double lm = log(m);
    return splitcomplex1_make((lp + lm) / 2, (lp - lm) / 2);
}

static inline splitcomplex1
splitcomplex1_pow(splitcomplex1 a, splitcomplex1 b)
{
    double ap = a.data[0] + a.data[1];
    double am = a.data[0] - a.data[1];
    double bp = b.data[0] + b.data[1];
    double bm = b.data[0] - b.data[1];
    double rp = pow(ap, bp);
    double rm = pow(am, bm);
    return splitcomplex1_make((rp + rm) / 2, (rp - rm) / 2);
}

static inline splitcomplex1
splitcomplex1_cos(splitcomplex1 a)
{
    return splitcomplex1_make(cos(a.data[0]) * cosh(a.data[1]),
                              -sin(a.data[0]) * sinh(a.data[1]));
}

static inline splitcomplex1
splitcomplex1_sin(splitcomplex1 a)
{
    return splitcomplex1_make(sin(a.data[0]) * cosh(a.data[1]),
                              cos(a.data[0]) * sinh(a.data[1]));
}

static inline splitcomplex1
splitcomplex1_tan(splitcomplex1 a)
{
    return splitcomplex1_div(splitcomplex1_sin(a), splitcomplex1_cos(a));
}

static inline splitcomplex1
splitcomplex1_cosh(splitcomplex1 a)
{
    return splitcomplex1_make(cosh(a.data[0]) * cosh(a.data[1]),
                              sinh(a.data[0]) * sinh(a.data[1]));
}

static inline splitcomplex1
splitcomplex1_sinh(splitcomplex1 a)
{
    return splitcomplex1_make(sinh(a.data[0]) * cosh(a.data[1]),
                              cosh(a.data[0]) * sinh(a.data[1]));
}

static inline splitcomplex1
splitcomplex1_tanh(splitcomplex1 a)
{
    return splitcomplex1_div(splitcomplex1_sinh(a), splitcomplex1_cosh(a));
}

static inline splitcomplex1
splitcomplex1_sqrt(splitcomplex1 a)
{
    double p = sqrt(a.data[0] + a.data[1]);
    double m = sqrt(a.data[0] - a.data[1]);
    return splitcomplex1_make((p + m) / 2, (p - m) / 2);
}


#endif /* ALGPY_SPLITCOMPLEX1_H */