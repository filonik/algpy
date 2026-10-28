/*
 * Arithmetic of split complex numbers (x + y*j), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_SPLITCOMPLEX1_H
#define ALGPY_SPLITCOMPLEX1_H

#include "../../fields/complex.h"

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

/* Evaluate the real components x+y and x-y independently.
 * Function domains are those of the real callback. */
static inline splitcomplex1
splitcomplex1_apply(splitcomplex1 a, double (*f)(double))
{
    double p = f(a.data[0] + a.data[1]);
    double q = f(a.data[0] - a.data[1]);
    return splitcomplex1_make((p + q)/2, (p - q)/2);
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
    int i, j;
    for (i = 0; i < ALGPY_SPLITCOMPLEX1_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX1_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
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
    return splitcomplex1_apply(a, exp);
}

static inline splitcomplex1
splitcomplex1_log(splitcomplex1 a)
{
    return splitcomplex1_apply(a, log);
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
    return splitcomplex1_apply(a, cos);
}

static inline splitcomplex1
splitcomplex1_sin(splitcomplex1 a)
{
    return splitcomplex1_apply(a, sin);
}

static inline splitcomplex1
splitcomplex1_tan(splitcomplex1 a)
{
    return splitcomplex1_apply(a, tan);
}

static inline splitcomplex1
splitcomplex1_cosh(splitcomplex1 a)
{
    return splitcomplex1_apply(a, cosh);
}

static inline splitcomplex1
splitcomplex1_sinh(splitcomplex1 a)
{
    return splitcomplex1_apply(a, sinh);
}

static inline splitcomplex1
splitcomplex1_tanh(splitcomplex1 a)
{
    return splitcomplex1_apply(a, tanh);
}

static inline splitcomplex1
splitcomplex1_sqrt(splitcomplex1 a)
{
    return splitcomplex1_apply(a, sqrt);
}

/* Common algebra API. */
static inline splitcomplex1
splitcomplex1_zero(void)
{
    return vec2_zero();
}

static inline int
splitcomplex1_is_real(splitcomplex1 a)
{
    return vec2_is_real(a);
}

static inline splitcomplex1
splitcomplex1_pos(splitcomplex1 a)
{
    return vec2_pos(a);
}

static inline splitcomplex1
splitcomplex1_neg(splitcomplex1 a)
{
    return vec2_neg(a);
}

static inline splitcomplex1
splitcomplex1_add(splitcomplex1 a, splitcomplex1 b)
{
    return vec2_add(a, b);
}

static inline splitcomplex1
splitcomplex1_sub(splitcomplex1 a, splitcomplex1 b)
{
    return vec2_sub(a, b);
}

static inline int
splitcomplex1_equal(splitcomplex1 a, splitcomplex1 b)
{
    return vec2_all_eq(a, b);
}

static inline int
splitcomplex1_not_equal(splitcomplex1 a, splitcomplex1 b)
{
    return !splitcomplex1_equal(a, b);
}

static inline splitcomplex1
splitcomplex1_acos(splitcomplex1 a)
{
    return splitcomplex1_apply(a, acos);
}

static inline splitcomplex1
splitcomplex1_asin(splitcomplex1 a)
{
    return splitcomplex1_apply(a, asin);
}

static inline splitcomplex1
splitcomplex1_atan(splitcomplex1 a)
{
    return splitcomplex1_apply(a, atan);
}

static inline splitcomplex1
splitcomplex1_acosh(splitcomplex1 a)
{
    return splitcomplex1_apply(a, acosh);
}

static inline splitcomplex1
splitcomplex1_asinh(splitcomplex1 a)
{
    return splitcomplex1_apply(a, asinh);
}

static inline splitcomplex1
splitcomplex1_atanh(splitcomplex1 a)
{
    return splitcomplex1_apply(a, atanh);
}

static inline splitcomplex1
splitcomplex1_square(splitcomplex1 a)
{
    return splitcomplex1_mul(a, a);
}

static inline splitcomplex1
splitcomplex1_radd(double a, splitcomplex1 b)
{
    return splitcomplex1_add(splitcomplex1_from_re(a), b);
}

static inline splitcomplex1
splitcomplex1_ladd(splitcomplex1 a, double b)
{
    return splitcomplex1_add(a, splitcomplex1_from_re(b));
}

static inline splitcomplex1
splitcomplex1_rsub(double a, splitcomplex1 b)
{
    return splitcomplex1_sub(splitcomplex1_from_re(a), b);
}

static inline splitcomplex1
splitcomplex1_lsub(splitcomplex1 a, double b)
{
    return splitcomplex1_sub(a, splitcomplex1_from_re(b));
}

static inline splitcomplex1
splitcomplex1_rmul(double a, splitcomplex1 b)
{
    return splitcomplex1_mul(splitcomplex1_from_re(a), b);
}

static inline splitcomplex1
splitcomplex1_lmul(splitcomplex1 a, double b)
{
    return splitcomplex1_mul(a, splitcomplex1_from_re(b));
}

static inline splitcomplex1
splitcomplex1_rdiv(double a, splitcomplex1 b)
{
    return splitcomplex1_div(splitcomplex1_from_re(a), b);
}

static inline splitcomplex1
splitcomplex1_ldiv(splitcomplex1 a, double b)
{
    return splitcomplex1_div(a, splitcomplex1_from_re(b));
}

static inline splitcomplex1
splitcomplex1_rpow(double a, splitcomplex1 b)
{
    return splitcomplex1_pow(splitcomplex1_from_re(a), b);
}

static inline splitcomplex1
splitcomplex1_lpow(splitcomplex1 a, double b)
{
    return splitcomplex1_pow(a, splitcomplex1_from_re(b));
}

static inline splitcomplex1
splitcomplex1_expm1(splitcomplex1 a)
{
    return splitcomplex1_apply(a, expm1);
}

static inline splitcomplex1
splitcomplex1_log1p(splitcomplex1 a)
{
    return splitcomplex1_apply(a, log1p);
}

static inline splitcomplex1
splitcomplex1_mix(splitcomplex1 a, splitcomplex1 b, double t)
{
    splitcomplex1 delta = splitcomplex1_mul(splitcomplex1_inv(a), b);
    return splitcomplex1_mul(a, splitcomplex1_lpow(delta, t));
}

static inline complex
splitcomplex1_norm(splitcomplex1 a)
{
    return complex_root_real(splitcomplex1_norm2(a), 2);
}

static inline void
splitcomplex1_minpoly(splitcomplex1 a, double result[3])
{
    result[0] = splitcomplex1_norm2(a);
    result[1] = -2*a.data[0];
    result[2] = 1;
}

static inline void
splitcomplex1_eivals(splitcomplex1 a, complex result[2])
{
    result[0] = complex_from_re(a.data[0] - a.data[1]);
    result[1] = complex_from_re(a.data[0] + a.data[1]);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
splitcomplex1_lmat(splitcomplex1 a, double result[ALGPY_SPLITCOMPLEX1_N*ALGPY_SPLITCOMPLEX1_N])
{
    splitcomplex1 coeffs[ALGPY_SPLITCOMPLEX1_N * ALGPY_SPLITCOMPLEX1_N];
    int row, col, i;
    splitcomplex1_coeffs(coeffs);
    for (row = 0; row < ALGPY_SPLITCOMPLEX1_N; row++) {
        for (col = 0; col < ALGPY_SPLITCOMPLEX1_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SPLITCOMPLEX1_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_SPLITCOMPLEX1_N + col].data[row];
            }
            result[row * ALGPY_SPLITCOMPLEX1_N + col] = value;
        }
    }
}

static inline void
splitcomplex1_rmat(splitcomplex1 a, double result[ALGPY_SPLITCOMPLEX1_N*ALGPY_SPLITCOMPLEX1_N])
{
    splitcomplex1 coeffs[ALGPY_SPLITCOMPLEX1_N * ALGPY_SPLITCOMPLEX1_N];
    int row, col, i;
    splitcomplex1_coeffs(coeffs);
    for (row = 0; row < ALGPY_SPLITCOMPLEX1_N; row++) {
        for (col = 0; col < ALGPY_SPLITCOMPLEX1_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SPLITCOMPLEX1_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_SPLITCOMPLEX1_N + i].data[col];
            }
            result[row * ALGPY_SPLITCOMPLEX1_N + col] = value;
        }
    }
}


#endif /* ALGPY_SPLITCOMPLEX1_H */
