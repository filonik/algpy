/*
 * Shared arithmetic of bi-complex numbers w + x*i0 + y*i1 + z*i0*i1
 * with i0*i0 = i1*i1 = -1 (so (i0*i1)^2 = +1), stored as a vec4.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX2_H
#define ALGPY_COMPLEX2_H

#include <math.h>

#include "../cayleydickson/complex.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_COMPLEX2_N ALGPY_VEC4_N

typedef vec4 complex2;

static inline complex2
complex2_make(double x, double x0, double x1, double x01)
{
    complex2 r = {{x, x0, x1, x01}};
    return r;
}

static inline complex2
complex2_from_re(double a)
{
    return complex2_make(a, 0, 0, 0);
}

static inline double
complex2_to_re(complex2 a)
{
    return a.data[0];
}

static inline complex2
complex2_one(void)
{
    return complex2_from_re(1);
}

static inline complex2
complex2_mul(complex2 x, complex2 y)
{
    const double *a = x.data, *b = y.data;
    return complex2_make(
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] + a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] - a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] - a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0]);
}

static inline void
complex2_basis(complex2 result[ALGPY_COMPLEX2_N])
{
    result[0] = complex2_make(1, 0, 0, 0);
    result[1] = complex2_make(0, 1, 0, 0);
    result[2] = complex2_make(0, 0, 1, 0);
    result[3] = complex2_make(0, 0, 0, 1);
}

static inline void
complex2_coeffs(
    complex2 result[ALGPY_COMPLEX2_N * ALGPY_COMPLEX2_N])
{
    complex2 basis[ALGPY_COMPLEX2_N];
    int i, j;
    complex2_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX2_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX2_N; j++) {
            result[i * ALGPY_COMPLEX2_N + j] =
                complex2_mul(basis[i], basis[j]);
        }
    }
}

/* Conjugation of the outer i1 part: (u + v*i1) -> (u - v*i1) */
static inline complex2
complex2_conjugate(complex2 a)
{
    return complex2_make(a.data[0], a.data[1], -a.data[2], -a.data[3]);
}

/* Norm a*conj(a) = u^2 + v^2 (a complex number) with a = u + v*i1. */
static inline complex
complex2_norm2(complex2 a)
{
    complex u = complex_make(a.data[0], a.data[1]);
    complex v = complex_make(a.data[2], a.data[3]);
    return vec2_add(complex_mul(u, u), complex_mul(v, v));
}

static inline double
complex2_norm4(complex2 a)
{
    return complex_norm2(complex2_norm2(a));
}

static inline complex
complex2_norm(complex2 a)
{
    return complex_root_real(complex2_norm4(a), 4);
}

static inline complex2
complex2_inv(complex2 a)
{
    complex ninv = complex_inv(complex2_norm2(a));
    complex u = complex_make(a.data[0], a.data[1]);
    complex v = complex_make(a.data[2], a.data[3]);
    complex ru = complex_mul(u, ninv);
    complex rv = complex_mul(vec2_neg(v), ninv);
    return complex2_make(ru.data[0], ru.data[1], rv.data[0], rv.data[1]);
}

static inline complex2
complex2_div(complex2 a, complex2 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return complex2_mul(a, complex2_inv(b));
}

/*
 * Apply a complex function f to a bi-complex number.
 *
 * Write a = u + v*i1 with u = w + x*i0 and v = y + z*i0 (u, v are complex
 * in i0). The algebra is commutative, so it splits into two copies of C,
 * one for each root of t^2 = -1 that i1 can map to:
 *   phi+ : i0 -> i, i1 -> +i,   a -> p = u + i*v = (w - z) + (x + y)*i
 *   phi- : i0 -> i, i1 -> -i,   a -> q = u - i*v = (w + z) + (x - y)*i
 * Then f(a) = phi+^-1(f(p)) and phi-^-1(f(q)) glued:
 *   u' = (f(p) + f(q)) / 2,  v' = (f(p) - f(q)) / (2i),
 * which gives the components computed below.
 */
static inline complex2
complex2_apply(complex2 a, complex (*f)(complex))
{
    complex p = f(complex_make(a.data[0] - a.data[3], a.data[1] + a.data[2]));
    complex q = f(complex_make(a.data[0] + a.data[3], a.data[1] - a.data[2]));
    double d0 = (q.data[0] - p.data[0]) / 2, d1 = (p.data[1] - q.data[1]) / 2;
    double d3 = (q.data[0] + p.data[0]) / 2, d2 = (p.data[1] + q.data[1]) / 2;
    return complex2_make(d3, d2, d1, d0);
}

static inline complex2
complex2_exp(complex2 a)
{
    return complex2_apply(a, complex_exp);
}

static inline complex2
complex2_log(complex2 a)
{
    return complex2_apply(a, complex_log);
}

static inline complex2
complex2_pow(complex2 a, complex2 b)
{
    return complex2_exp(complex2_mul(complex2_log(a), b));
}

static inline complex2
complex2_sqrt(complex2 a)
{
    return complex2_apply(a, complex_sqrt);
}

/* Minimal polynomial coefficients (degree 4, monic), lowest degree first. */
static inline void
complex2_minpoly(complex2 a, double result[5])
{
    double x0 = a.data[0], x1 = a.data[1], x2 = a.data[2], x3 = a.data[3];
    double x02 = x0*x0, x12 = x1*x1, x22 = x2*x2, x32 = x3*x3;
    result[4] = 1;
    result[3] = -4*x0;
    result[2] = -2*x32 + 2*x12 + 2*x22 + 6*x02;
    result[1] = 4*x0*x32 - 8*x1*x2*x3 - 4*x0*x12 - 4*x0*x22 - 4*x02*x0;
    result[0] = x32*x32 + 2*x12*x32 + x12*x12 + 2*x22*x32 - 2*x02*x32
        + 8*x0*x1*x2*x3 - 2*x12*x22 + 2*x02*x12 + x22*x22 + 2*x02*x22
        + x02*x02;
}

/* Eigenvalues (x0 -/+ x3) +/- i*(x1 +/- x2) as complex numbers. */
static inline void
complex2_eivals(complex2 a, complex result[4])
{
    double x0 = a.data[0], x1 = a.data[1], x2 = a.data[2], x3 = a.data[3];
    result[0] = complex_make(x0 - x3, -x1 - x2);
    result[1] = complex_make(x0 - x3, x1 + x2);
    result[2] = complex_make(x0 + x3, x1 - x2);
    result[3] = complex_make(x0 + x3, -x1 + x2);
}

#endif /* ALGPY_COMPLEX2_H */
