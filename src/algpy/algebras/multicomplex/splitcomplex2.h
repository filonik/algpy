/*
 * Arithmetic of bi-split-complex numbers x + x0*e0 + x1*e1 + x01*e0*e1 with
 * e0*e0 = e1*e1 = +1, stored as a vec4.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_SPLITCOMPLEX2_H
#define ALGPY_SPLITCOMPLEX2_H

#include <math.h>

#include "../../fields/complex.h"
#include "splitcomplex1.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_SPLITCOMPLEX2_N ALGPY_VEC4_N

typedef vec4 splitcomplex2;

static inline splitcomplex2
splitcomplex2_make(double x, double x0, double x1, double x01)
{
    splitcomplex2 r = {{x, x0, x1, x01}};
    return r;
}

static inline splitcomplex2
splitcomplex2_from_re(double a)
{
    return splitcomplex2_make(a, 0, 0, 0);
}

static inline double
splitcomplex2_to_re(splitcomplex2 a)
{
    return a.data[0];
}

static inline splitcomplex2
splitcomplex2_one(void)
{
    return splitcomplex2_from_re(1);
}

static inline splitcomplex2
splitcomplex2_mul(splitcomplex2 x, splitcomplex2 y)
{
    const double *a = x.data, *b = y.data;
    return splitcomplex2_make(
        a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] + a[3]*b[2],
        a[0]*b[2] + a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0]);
}

static inline void
splitcomplex2_basis(splitcomplex2 result[ALGPY_SPLITCOMPLEX2_N])
{
    int i, j;
    for (i = 0; i < ALGPY_SPLITCOMPLEX2_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX2_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
splitcomplex2_coeffs(
    splitcomplex2 result[ALGPY_SPLITCOMPLEX2_N * ALGPY_SPLITCOMPLEX2_N])
{
    splitcomplex2 basis[ALGPY_SPLITCOMPLEX2_N];
    int i, j;
    splitcomplex2_basis(basis);
    for (i = 0; i < ALGPY_SPLITCOMPLEX2_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX2_N; j++) {
            result[i * ALGPY_SPLITCOMPLEX2_N + j] =
                splitcomplex2_mul(basis[i], basis[j]);
        }
    }
}

/* Conjugation of the outer i1 part: (u + v*i1) -> (u - v*i1) */
static inline splitcomplex2
splitcomplex2_conjugate(splitcomplex2 a)
{
    return splitcomplex2_make(a.data[0], a.data[1], -a.data[2], -a.data[3]);
}

/* Norm a*conj(a) = u^2 - v^2 (a split-complex number) with a = u + v*i1. */
static inline splitcomplex1
splitcomplex2_norm2(splitcomplex2 a)
{
    splitcomplex1 u = splitcomplex1_make(a.data[0], a.data[1]);
    splitcomplex1 v = splitcomplex1_make(a.data[2], a.data[3]);
    return vec2_sub(splitcomplex1_mul(u, u), splitcomplex1_mul(v, v));
}

static inline double
splitcomplex2_norm4(splitcomplex2 a)
{
    return splitcomplex1_norm2(splitcomplex2_norm2(a));
}

static inline complex
splitcomplex2_norm(splitcomplex2 a)
{
    return complex_root_real(splitcomplex2_norm4(a), 4);
}

static inline splitcomplex2
splitcomplex2_inv(splitcomplex2 a)
{
    splitcomplex1 ninv = splitcomplex1_inv(splitcomplex2_norm2(a));
    splitcomplex1 u = splitcomplex1_make(a.data[0], a.data[1]);
    splitcomplex1 v = splitcomplex1_make(a.data[2], a.data[3]);
    splitcomplex1 ru = splitcomplex1_mul(u, ninv);
    splitcomplex1 rv = splitcomplex1_mul(vec2_neg(v), ninv);
    return splitcomplex2_make(ru.data[0], ru.data[1], rv.data[0], rv.data[1]);
}

static inline splitcomplex2
splitcomplex2_div(splitcomplex2 a, splitcomplex2 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return splitcomplex2_mul(a, splitcomplex2_inv(b));
}

/*
 * The algebra is isomorphic to R^4 through the four idempotent components
 * l[2*i + j] = x + s_i*x0 + s_j*x1 + s_i*s_j*x01 with s_0 = +1, s_1 = -1.
 * Multiplication is componentwise there, so smooth real functions are applied
 * to each component.
 */
static inline void
splitcomplex2_components(splitcomplex2 a, double l[4])
{
    const double *x = a.data;
    l[0] = x[0] + x[1] + x[2] + x[3];
    l[1] = x[0] + x[1] - x[2] - x[3];
    l[2] = x[0] - x[1] + x[2] - x[3];
    l[3] = x[0] - x[1] - x[2] + x[3];
}

static inline splitcomplex2
splitcomplex2_from_components(const double l[4])
{
    return splitcomplex2_make(
        (l[0] + l[1] + l[2] + l[3]) / 4,
        (l[0] + l[1] - l[2] - l[3]) / 4,
        (l[0] - l[1] + l[2] - l[3]) / 4,
        (l[0] - l[1] - l[2] + l[3]) / 4);
}

/* For a = u + v*e with e*e = +1, apply f to u+v and u-v. */
static inline splitcomplex2
splitcomplex2_apply(splitcomplex2 a, splitcomplex1 (*f)(splitcomplex1))
{
    splitcomplex1 p = splitcomplex1_make(
        a.data[0] + a.data[2],
        a.data[1] + a.data[3]);
    splitcomplex1 q = splitcomplex1_make(
        a.data[0] - a.data[2],
        a.data[1] - a.data[3]);
    p = f(p);
    q = f(q);
    return splitcomplex2_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2);
}

static inline splitcomplex2
splitcomplex2_exp(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_exp);
}

static inline splitcomplex2
splitcomplex2_log(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_log);
}

static inline splitcomplex2
splitcomplex2_pow(splitcomplex2 a, splitcomplex2 b)
{
    splitcomplex1 ap = splitcomplex1_make(
        a.data[0] + a.data[2],
        a.data[1] + a.data[3]);
    splitcomplex1 aq = splitcomplex1_make(
        a.data[0] - a.data[2],
        a.data[1] - a.data[3]);
    splitcomplex1 bp = splitcomplex1_make(
        b.data[0] + b.data[2],
        b.data[1] + b.data[3]);
    splitcomplex1 bq = splitcomplex1_make(
        b.data[0] - b.data[2],
        b.data[1] - b.data[3]);
    splitcomplex1 p = splitcomplex1_pow(ap, bp);
    splitcomplex1 q = splitcomplex1_pow(aq, bq);
    return splitcomplex2_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2);
}

static inline splitcomplex2
splitcomplex2_sqrt(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_sqrt);
}

static inline splitcomplex2
splitcomplex2_cos(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_cos);
}

static inline splitcomplex2
splitcomplex2_sin(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_sin);
}

static inline splitcomplex2
splitcomplex2_tan(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_tan);
}

static inline splitcomplex2
splitcomplex2_cosh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_cosh);
}

static inline splitcomplex2
splitcomplex2_sinh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_sinh);
}

static inline splitcomplex2
splitcomplex2_tanh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_tanh);
}

static inline splitcomplex2
splitcomplex2_radd(double a, splitcomplex2 b) {
  return vec4_add(splitcomplex2_from_re(a), b);
}

static inline splitcomplex2
splitcomplex2_ladd(splitcomplex2 a, double b) {
  return vec4_add(a, splitcomplex2_from_re(b));
}

static inline splitcomplex2
splitcomplex2_rsub(double a, splitcomplex2 b) {
  return vec4_sub(splitcomplex2_from_re(a), b);
}

static inline splitcomplex2
splitcomplex2_lsub(splitcomplex2 a, double b) {
  return vec4_sub(a, splitcomplex2_from_re(b));
}

static inline splitcomplex2
splitcomplex2_rmul(double a, splitcomplex2 b) {
  return splitcomplex2_mul(splitcomplex2_from_re(a), b);
}

static inline splitcomplex2
splitcomplex2_lmul(splitcomplex2 a, double b) {
  return splitcomplex2_mul(a, splitcomplex2_from_re(b));
}

static inline splitcomplex2
splitcomplex2_rdiv(double a, splitcomplex2 b) {
  return splitcomplex2_div(splitcomplex2_from_re(a), b);
}

static inline splitcomplex2
splitcomplex2_ldiv(splitcomplex2 a, double b) {
  return splitcomplex2_div(a, splitcomplex2_from_re(b));
}

static inline splitcomplex2
splitcomplex2_rpow(double a, splitcomplex2 b) {
  return splitcomplex2_pow(splitcomplex2_from_re(a), b);
}

static inline splitcomplex2
splitcomplex2_lpow(splitcomplex2 a, double b) {
  return splitcomplex2_pow(a, splitcomplex2_from_re(b));
}

static inline splitcomplex2
splitcomplex2_mix(splitcomplex2 a, splitcomplex2 b, double t)
{
    splitcomplex2 delta = splitcomplex2_mul(splitcomplex2_inv(a), b);
    return splitcomplex2_mul(a, splitcomplex2_lpow(delta, t));
}

static inline void
splitcomplex2_minpoly(splitcomplex2 a, double result[5])
{
    double l[4];
    int i, j;
    splitcomplex2_components(a, l);
    result[0] = 1;
    for (j = 1; j < 5; j++) {
        result[j] = 0;
    }
    for (i = 0; i < 4; i++) {
        for (j = i + 1; j > 0; j--) {
            result[j] = result[j - 1] - l[i] * result[j];
        }
        result[0] = -l[i] * result[0];
    }
}

static inline void
splitcomplex2_eivals(splitcomplex2 a, complex result[4])
{
    double l[4];
    int i;
    splitcomplex2_components(a, l);
    for (i = 0; i < 4; i++) {
        result[i] = complex_make(l[i], 0);
    }
}

/* Common algebra API. */
static inline splitcomplex2
splitcomplex2_zero(void)
{
    return vec4_zero();
}

static inline int
splitcomplex2_is_real(splitcomplex2 a)
{
    return vec4_is_real(a);
}

static inline splitcomplex2
splitcomplex2_pos(splitcomplex2 a)
{
    return vec4_pos(a);
}

static inline splitcomplex2
splitcomplex2_neg(splitcomplex2 a)
{
    return vec4_neg(a);
}

static inline splitcomplex2
splitcomplex2_add(splitcomplex2 a, splitcomplex2 b)
{
    return vec4_add(a, b);
}

static inline splitcomplex2
splitcomplex2_sub(splitcomplex2 a, splitcomplex2 b)
{
    return vec4_sub(a, b);
}

static inline int
splitcomplex2_equal(splitcomplex2 a, splitcomplex2 b)
{
    return vec4_all_eq(a, b);
}

static inline int
splitcomplex2_not_equal(splitcomplex2 a, splitcomplex2 b)
{
    return !splitcomplex2_equal(a, b);
}

static inline splitcomplex2
splitcomplex2_re(splitcomplex2 a)
{
    return splitcomplex2_from_re(a.data[0]);
}

static inline splitcomplex2
splitcomplex2_im(splitcomplex2 a)
{
    return splitcomplex2_sub(a, splitcomplex2_re(a));
}

static inline splitcomplex2
splitcomplex2_acos(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_acos);
}

static inline splitcomplex2
splitcomplex2_asin(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_asin);
}

static inline splitcomplex2
splitcomplex2_atan(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_atan);
}

static inline splitcomplex2
splitcomplex2_acosh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_acosh);
}

static inline splitcomplex2
splitcomplex2_asinh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_asinh);
}

static inline splitcomplex2
splitcomplex2_atanh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_atanh);
}

static inline splitcomplex2
splitcomplex2_square(splitcomplex2 a)
{
    return splitcomplex2_mul(a, a);
}

static inline splitcomplex2
splitcomplex2_expm1(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_expm1);
}

static inline splitcomplex2
splitcomplex2_log1p(splitcomplex2 a)
{
    return splitcomplex2_apply(a, splitcomplex1_log1p);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
splitcomplex2_lmat(splitcomplex2 a, double result[ALGPY_SPLITCOMPLEX2_N*ALGPY_SPLITCOMPLEX2_N])
{
    splitcomplex2 coeffs[ALGPY_SPLITCOMPLEX2_N * ALGPY_SPLITCOMPLEX2_N];
    int row, col, i;
    splitcomplex2_coeffs(coeffs);
    for (row = 0; row < ALGPY_SPLITCOMPLEX2_N; row++) {
        for (col = 0; col < ALGPY_SPLITCOMPLEX2_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SPLITCOMPLEX2_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_SPLITCOMPLEX2_N + col].data[row];
            }
            result[row * ALGPY_SPLITCOMPLEX2_N + col] = value;
        }
    }
}

static inline void
splitcomplex2_rmat(splitcomplex2 a, double result[ALGPY_SPLITCOMPLEX2_N*ALGPY_SPLITCOMPLEX2_N])
{
    splitcomplex2 coeffs[ALGPY_SPLITCOMPLEX2_N * ALGPY_SPLITCOMPLEX2_N];
    int row, col, i;
    splitcomplex2_coeffs(coeffs);
    for (row = 0; row < ALGPY_SPLITCOMPLEX2_N; row++) {
        for (col = 0; col < ALGPY_SPLITCOMPLEX2_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SPLITCOMPLEX2_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_SPLITCOMPLEX2_N + i].data[col];
            }
            result[row * ALGPY_SPLITCOMPLEX2_N + col] = value;
        }
    }
}


#endif /* ALGPY_SPLITCOMPLEX2_H */
