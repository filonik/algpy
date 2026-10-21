/*
 * Core arithmetic of the multicomplex numbers x + y*e with e*e = k.
 *
 * Pure C with no dependency on Python or NumPy. The parameter k is passed
 * explicitly to every operation that depends on it.
 */
#ifndef ALGPY_MULTICOMPLEX1_H
#define ALGPY_MULTICOMPLEX1_H

#define ALGPY_MULTICOMPLEX1_N 2

#include <math.h>

#include "../../complex/complex1.h"
#include "../../../vectorspaces/vec1.h"
#include "../../../vectorspaces/vec2.h"

typedef vec2 multicomplex1;

static inline multicomplex1
multicomplex1_make(double x, double y)
{
    multicomplex1 r = {{x, y}};
    return r;
}

static inline multicomplex1
multicomplex1_re(multicomplex1 a)
{
    return multicomplex1_make(a.data[0], 0);
}

static inline multicomplex1
multicomplex1_im(multicomplex1 a)
{
    return multicomplex1_make(0, a.data[1]);
}

static inline multicomplex1
multicomplex1_from_re(double a)
{
    return multicomplex1_make(a, 0);
}

static inline double
multicomplex1_to_re(multicomplex1 a)
{
    return a.data[0];
}

static inline multicomplex1
multicomplex1_one(void)
{
    return multicomplex1_from_re(1);
}

static inline multicomplex1
multicomplex1_mul(vec1 k, multicomplex1 x, multicomplex1 y)
{
    if (vec2_is_real(y)) {
        return vec2_scale(x, y.data[0]);
    }
    if (vec2_is_real(x)) {
        return vec2_scale(y, x.data[0]);
    }
    return multicomplex1_make(
        x.data[0]*y.data[0] + k.data[0]*x.data[1]*y.data[1],
        x.data[0]*y.data[1] + x.data[1]*y.data[0]);
}

static inline multicomplex1
multicomplex1_conjugate(multicomplex1 a)
{
    return multicomplex1_make(a.data[0], -a.data[1]);
}

/* Norm a*conj(a) = x^2 - k*y^2 */
static inline double
multicomplex1_norm2(vec1 k, multicomplex1 a)
{
    return a.data[0]*a.data[0] - k.data[0]*a.data[1]*a.data[1];
}

static inline multicomplex1
multicomplex1_inv(vec1 k, multicomplex1 a)
{
    double n = multicomplex1_norm2(k, a);
    return multicomplex1_make(a.data[0] / n, -a.data[1] / n);
}

static inline multicomplex1
multicomplex1_div(vec1 k, multicomplex1 a, multicomplex1 b)
{
    if (vec2_is_real(b)) {
        return vec2_div_scalar(a, b.data[0]);
    }
    return multicomplex1_mul(k, a, multicomplex1_inv(k, b));
}

static inline void
multicomplex1_basis(multicomplex1 result[ALGPY_MULTICOMPLEX1_N])
{
    result[0] = multicomplex1_make(1, 0);
    result[1] = multicomplex1_make(0, 1);
}

static inline void
multicomplex1_coeffs(
    vec1 k,
    multicomplex1 result[ALGPY_MULTICOMPLEX1_N * ALGPY_MULTICOMPLEX1_N])
{
    multicomplex1 basis[ALGPY_MULTICOMPLEX1_N];
    int i, j;
    multicomplex1_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX1_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX1_N; j++) {
            result[i * ALGPY_MULTICOMPLEX1_N + j] =
                multicomplex1_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
multicomplex1_minpoly(vec1 k, multicomplex1 a, double result[3])
{
    result[0] = multicomplex1_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/* Eigenvalues x -/+ y*sqrt(k), as complex numbers. */
static inline void
multicomplex1_eivals(vec1 k, multicomplex1 a, complex1 result[2])
{
    complex1 s = complex1_sqrt_real(k.data[0]);
    complex1 t = complex1_make(a.data[1], 0);
    complex1 d = complex1_mul(t, s);
    result[0] = complex1_make(a.data[0] - d.data[0], -d.data[1]);
    result[1] = complex1_make(a.data[0] + d.data[0], +d.data[1]);
}

#endif /* ALGPY_MULTICOMPLEX1_H */
