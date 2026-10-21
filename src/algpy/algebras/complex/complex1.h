/*
 * Shared arithmetic of complex numbers (x + y*i), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX1_H
#define ALGPY_COMPLEX1_H

#define ALGPY_COMPLEX1_N 2

#include <math.h>

#include "../../vectorspaces/vec2.h"

typedef vec2 complex1;

static inline complex1
complex1_make(double x, double y)
{
    complex1 r = {{x, y}};
    return r;
}

static inline complex1
complex1_re(complex1 a)
{
    return complex1_make(a.data[0], 0);
}

static inline complex1
complex1_im(complex1 a)
{
    return complex1_make(0, a.data[1]);
}

static inline complex1
complex1_from_re(double a)
{
    return complex1_make(a, 0);
}

static inline double
complex1_to_re(complex1 a)
{
    return a.data[0];
}

static inline complex1
complex1_one(void)
{
    return complex1_from_re(1);
}

static inline complex1
complex1_mul(complex1 x, complex1 y)
{
    return complex1_make(
        x.data[0]*y.data[0] - x.data[1]*y.data[1],
        x.data[0]*y.data[1] + x.data[1]*y.data[0]);
}

static inline void
complex1_basis(complex1 result[ALGPY_COMPLEX1_N])
{
    result[0] = complex1_make(1, 0);
    result[1] = complex1_make(0, 1);
}

static inline void
complex1_coeffs(
    complex1 result[ALGPY_COMPLEX1_N * ALGPY_COMPLEX1_N])
{
    complex1 basis[ALGPY_COMPLEX1_N];
    int i, j;
    complex1_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX1_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX1_N; j++) {
            result[i * ALGPY_COMPLEX1_N + j] =
                complex1_mul(basis[i], basis[j]);
        }
    }
}

static inline complex1
complex1_conjugate(complex1 a)
{
    return complex1_make(a.data[0], -a.data[1]);
}

/* Norm a*conj(a) = x^2 + y^2 */
static inline double
complex1_norm2(complex1 a)
{
    return a.data[0]*a.data[0] + a.data[1]*a.data[1];
}

static inline complex1
complex1_inv(complex1 a)
{
    double n = complex1_norm2(a);
    return complex1_make(a.data[0] / n, -a.data[1] / n);
}

static inline complex1
complex1_div(complex1 a, complex1 b)
{
    return complex1_mul(a, complex1_inv(b));
}

static inline complex1
complex1_exp(complex1 a)
{
    double e = exp(a.data[0]);
    return complex1_make(e * cos(a.data[1]), e * sin(a.data[1]));
}

static inline complex1
complex1_log(complex1 a)
{
    double r = hypot(a.data[0], a.data[1]);
    return complex1_make(r == 0 ? -INFINITY : log(r),
                         atan2(a.data[1], a.data[0]));
}

static inline complex1
complex1_pow(complex1 a, complex1 b)
{
    return complex1_exp(complex1_mul(complex1_log(a), b));
}

static inline complex1
complex1_cos(complex1 a)
{
    return complex1_make(cos(a.data[0]) * cosh(a.data[1]),
                         -sin(a.data[0]) * sinh(a.data[1]));
}

static inline complex1
complex1_sin(complex1 a)
{
    return complex1_make(sin(a.data[0]) * cosh(a.data[1]),
                         cos(a.data[0]) * sinh(a.data[1]));
}

static inline complex1
complex1_tan(complex1 a)
{
    return complex1_div(complex1_sin(a), complex1_cos(a));
}

static inline complex1
complex1_cosh(complex1 a)
{
    return complex1_make(cosh(a.data[0]) * cos(a.data[1]),
                         sinh(a.data[0]) * sin(a.data[1]));
}

static inline complex1
complex1_sinh(complex1 a)
{
    return complex1_make(sinh(a.data[0]) * cos(a.data[1]),
                         cosh(a.data[0]) * sin(a.data[1]));
}

static inline complex1
complex1_tanh(complex1 a)
{
    return complex1_div(complex1_sinh(a), complex1_cosh(a));
}

static inline complex1
complex1_sqrt(complex1 a)
{
    double r = hypot(a.data[0], a.data[1]);
    double re = sqrt(fmax(0, (r + a.data[0]) / 2));
    double im = sqrt(fmax(0, (r - a.data[0]) / 2));
    return complex1_make(re, copysign(im, a.data[1]));
}

/* Square root of a real number: sqrt(k) or i*sqrt(-k). */
static inline complex1
complex1_sqrt_real(double k)
{
    if (k < 0) {
        return complex1_make(0, sqrt(-k));
    }
    return complex1_make(sqrt(k), 0);
}


#endif /* ALGPY_COMPLEX1_H */