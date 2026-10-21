/*
 * Arithmetic of dual numbers (x + y*e), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_DUAL1_H
#define ALGPY_DUAL1_H

#define ALGPY_DUAL1_N 2

#include <math.h>

#include "../../vectorspaces/vec2.h"

typedef vec2 dual1;

static inline dual1
dual1_make(double x, double y)
{
    dual1 r = {{x, y}};
    return r;
}

static inline dual1
dual1_re(dual1 a)
{
    return dual1_make(a.data[0], 0);
}

static inline dual1
dual1_im(dual1 a)
{
    return dual1_make(0, a.data[1]);
}

static inline dual1
dual1_from_re(double a)
{
    return dual1_make(a, 0);
}

static inline double
dual1_to_re(dual1 a)
{
    return a.data[0];
}

static inline dual1
dual1_one(void)
{
    return dual1_from_re(1);
}

static inline dual1
dual1_mul(dual1 x, dual1 y)
{
    return dual1_make(
        x.data[0]*y.data[0],
        x.data[0]*y.data[1] + x.data[1]*y.data[0]);
}

static inline void
dual1_basis(dual1 result[ALGPY_DUAL1_N])
{
    result[0] = dual1_make(1, 0);
    result[1] = dual1_make(0, 1);
}

static inline void
dual1_coeffs(
    dual1 result[ALGPY_DUAL1_N * ALGPY_DUAL1_N])
{
    dual1 basis[ALGPY_DUAL1_N];
    int i, j;
    dual1_basis(basis);
    for (i = 0; i < ALGPY_DUAL1_N; i++) {
        for (j = 0; j < ALGPY_DUAL1_N; j++) {
            result[i * ALGPY_DUAL1_N + j] =
                dual1_mul(basis[i], basis[j]);
        }
    }
}

static inline dual1
dual1_conjugate(dual1 a)
{
    return dual1_make(a.data[0], -a.data[1]);
}

/* Norm a*conj(a) = x^2 */
static inline double
dual1_norm2(dual1 a)
{
    return a.data[0]*a.data[0];
}

static inline dual1
dual1_inv(dual1 a)
{
    double n = dual1_norm2(a);
    return dual1_make(a.data[0] / n, -a.data[1] / n);
}

static inline dual1
dual1_div(dual1 a, dual1 b)
{
    return dual1_mul(a, dual1_inv(b));
}

static inline dual1
dual1_exp(dual1 a)
{
    double e = exp(a.data[0]);
    return dual1_make(e, e * a.data[1]);
}

static inline dual1
dual1_log(dual1 a)
{
    return dual1_make(log(a.data[0]), a.data[1] / a.data[0]);
}

static inline dual1
dual1_pow(dual1 a, dual1 b)
{
    return dual1_exp(dual1_mul(dual1_log(a), b));
}

static inline dual1
dual1_cos(dual1 a)
{
    return dual1_make(cos(a.data[0]), -sin(a.data[0]) * a.data[1]);
}

static inline dual1
dual1_sin(dual1 a)
{
    return dual1_make(sin(a.data[0]), cos(a.data[0]) * a.data[1]);
}

static inline dual1
dual1_tan(dual1 a)
{
    return dual1_div(dual1_sin(a), dual1_cos(a));
}

static inline dual1
dual1_cosh(dual1 a)
{
    return dual1_make(cosh(a.data[0]), sinh(a.data[0]) * a.data[1]);
}

static inline dual1
dual1_sinh(dual1 a)
{
    return dual1_make(sinh(a.data[0]), cosh(a.data[0]) * a.data[1]);
}

static inline dual1
dual1_tanh(dual1 a)
{
    return dual1_div(dual1_sinh(a), dual1_cosh(a));
}

static inline dual1
dual1_sqrt(dual1 a)
{
    double r = sqrt(a.data[0]);
    return dual1_make(r, a.data[1] / (2 * r));
}


#endif /* ALGPY_DUAL1_H */