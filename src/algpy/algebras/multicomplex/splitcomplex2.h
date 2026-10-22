/*
 * Arithmetic of bi-split-complex numbers x + x0*e0 + x1*e1 + x01*e0*e1 with
 * e0*e0 = e1*e1 = +1, stored as a vec4.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_SPLITCOMPLEX2_H
#define ALGPY_SPLITCOMPLEX2_H

#include <math.h>

#include "../cayleydickson/complex.h"
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
    result[0] = splitcomplex2_make(1, 0, 0, 0);
    result[1] = splitcomplex2_make(0, 1, 0, 0);
    result[2] = splitcomplex2_make(0, 0, 1, 0);
    result[3] = splitcomplex2_make(0, 0, 0, 1);
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

static inline splitcomplex2
splitcomplex2_apply(splitcomplex2 a, double (*f)(double))
{
    double l[4];
    int i;
    splitcomplex2_components(a, l);
    for (i = 0; i < 4; i++) {
        l[i] = f(l[i]);
    }
    return splitcomplex2_from_components(l);
}

static inline splitcomplex2
splitcomplex2_exp(splitcomplex2 a)
{
    return splitcomplex2_apply(a, exp);
}

static inline splitcomplex2
splitcomplex2_log(splitcomplex2 a)
{
    return splitcomplex2_apply(a, log);
}

static inline splitcomplex2
splitcomplex2_pow(splitcomplex2 a, splitcomplex2 b)
{
    return splitcomplex2_exp(splitcomplex2_mul(splitcomplex2_log(a), b));
}

static inline splitcomplex2
splitcomplex2_sqrt(splitcomplex2 a)
{
    return splitcomplex2_apply(a, sqrt);
}

static inline splitcomplex2
splitcomplex2_cos(splitcomplex2 a)
{
    return splitcomplex2_apply(a, cos);
}

static inline splitcomplex2
splitcomplex2_sin(splitcomplex2 a)
{
    return splitcomplex2_apply(a, sin);
}

static inline splitcomplex2
splitcomplex2_tan(splitcomplex2 a)
{
    return splitcomplex2_apply(a, tan);
}

static inline splitcomplex2
splitcomplex2_cosh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, cosh);
}

static inline splitcomplex2
splitcomplex2_sinh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, sinh);
}

static inline splitcomplex2
splitcomplex2_tanh(splitcomplex2 a)
{
    return splitcomplex2_apply(a, tanh);
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

#endif /* ALGPY_SPLITCOMPLEX2_H */
