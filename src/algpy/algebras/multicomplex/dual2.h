/*
 * Arithmetic of bi-dual numbers x + x0*e0 + x1*e1 + x01*e0*e1 with
 * e0*e0 = e1*e1 = 0, stored as a vec4.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_DUAL2_H
#define ALGPY_DUAL2_H

#include <math.h>

#include "../cayleydickson/complex.h"
#include "dual1.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_DUAL2_N ALGPY_VEC4_N

typedef vec4 dual2;

static inline dual2
dual2_make(double x, double x0, double x1, double x01)
{
    dual2 r = {{x, x0, x1, x01}};
    return r;
}

static inline dual2
dual2_from_re(double a)
{
    return dual2_make(a, 0, 0, 0);
}

static inline double
dual2_to_re(dual2 a)
{
    return a.data[0];
}

static inline dual2
dual2_one(void)
{
    return dual2_from_re(1);
}

static inline dual2
dual2_mul(dual2 x, dual2 y)
{
    const double *a = x.data, *b = y.data;
    return dual2_make(
        a[0]*b[0],
        a[0]*b[1] + a[1]*b[0],
        a[0]*b[2] + a[2]*b[0],
        a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0]);
}

static inline void
dual2_basis(dual2 result[ALGPY_DUAL2_N])
{
    result[0] = dual2_make(1, 0, 0, 0);
    result[1] = dual2_make(0, 1, 0, 0);
    result[2] = dual2_make(0, 0, 1, 0);
    result[3] = dual2_make(0, 0, 0, 1);
}

static inline void
dual2_coeffs(
    dual2 result[ALGPY_DUAL2_N * ALGPY_DUAL2_N])
{
    dual2 basis[ALGPY_DUAL2_N];
    int i, j;
    dual2_basis(basis);
    for (i = 0; i < ALGPY_DUAL2_N; i++) {
        for (j = 0; j < ALGPY_DUAL2_N; j++) {
            result[i * ALGPY_DUAL2_N + j] =
                dual2_mul(basis[i], basis[j]);
        }
    }
}

/* Conjugation of the outer i1 part: (u + v*i1) -> (u - v*i1) */
static inline dual2
dual2_conjugate(dual2 a)
{
    return dual2_make(a.data[0], a.data[1], -a.data[2], -a.data[3]);
}

/* Norm a*conj(a) = u^2 (a dual number) with a = u + v*i1, i1^2 = 0. */
static inline dual1
dual2_norm2(dual2 a)
{
    dual1 u = dual1_make(a.data[0], a.data[1]);
    return dual1_mul(u, u);
}

static inline double
dual2_norm4(dual2 a)
{
    return dual1_norm2(dual2_norm2(a));
}

static inline complex
dual2_norm(dual2 a)
{
    return complex_root_real(dual2_norm4(a), 4);
}

static inline dual2
dual2_inv(dual2 a)
{
    dual1 ninv = dual1_inv(dual2_norm2(a));
    dual1 u = dual1_make(a.data[0], a.data[1]);
    dual1 v = dual1_make(a.data[2], a.data[3]);
    dual1 ru = dual1_mul(u, ninv);
    dual1 rv = dual1_mul(vec2_neg(v), ninv);
    return dual2_make(ru.data[0], ru.data[1], rv.data[0], rv.data[1]);
}

static inline dual2
dual2_div(dual2 a, dual2 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return dual2_mul(a, dual2_inv(b));
}

/*
 * Apply a smooth real function to a number x + x0*e0 + x1*e1 + x01*e0*e1
 * (e0^2 = e1^2 = 0) given f(x), f'(x) and f''(x):
 *   f(a) = f + f'*(x0*e0 + x1*e1 + x01*e0*e1) + f''*x0*x1*e0*e1.
 */
static inline dual2
dual2_apply(dual2 a, double f0, double f1, double f2)
{
    return dual2_make(
        f0,
        f1*a.data[1],
        f1*a.data[2],
        f1*a.data[3] + f2*a.data[1]*a.data[2]);
}

static inline dual2
dual2_exp(dual2 a)
{
    double e = exp(a.data[0]);
    return dual2_apply(a, e, e, e);
}

static inline dual2
dual2_log(dual2 a)
{
    double x = a.data[0];
    return dual2_apply(a, log(x), 1 / x, -1 / (x*x));
}

static inline dual2
dual2_pow(dual2 a, dual2 b)
{
    return dual2_exp(dual2_mul(dual2_log(a), b));
}

static inline dual2
dual2_sqrt(dual2 a)
{
    double r = sqrt(a.data[0]);
    return dual2_apply(a, r, 1 / (2*r), -1 / (4*r*a.data[0]));
}

static inline dual2
dual2_radd(double a, dual2 b) {
  return vec4_add(dual2_from_re(a), b);
}

static inline dual2
dual2_ladd(dual2 a, double b) {
  return vec4_add(a, dual2_from_re(b));
}

static inline dual2
dual2_rsub(double a, dual2 b) {
  return vec4_sub(dual2_from_re(a), b);
}

static inline dual2
dual2_lsub(dual2 a, double b) {
  return vec4_sub(a, dual2_from_re(b));
}

static inline dual2
dual2_rmul(double a, dual2 b) {
  return dual2_mul(dual2_from_re(a), b);
}

static inline dual2
dual2_lmul(dual2 a, double b) {
  return dual2_mul(a, dual2_from_re(b));
}

static inline dual2
dual2_rdiv(double a, dual2 b) {
  return dual2_div(dual2_from_re(a), b);
}

static inline dual2
dual2_ldiv(dual2 a, double b) {
  return dual2_div(a, dual2_from_re(b));
}

static inline dual2
dual2_rpow(double a, dual2 b) {
  return dual2_pow(dual2_from_re(a), b);
}

static inline dual2
dual2_lpow(dual2 a, double b) {
  return dual2_pow(a, dual2_from_re(b));
}

static inline dual2
dual2_mix(dual2 a, dual2 b, double t)
{
    dual2 delta = dual2_mul(dual2_inv(a), b);
    return dual2_mul(a, dual2_lpow(delta, t));
}

static inline void
dual2_minpoly(dual2 a, double result[5])
{
    double x0 = a.data[0];
    result[4] = 1;
    result[3] = -4*x0;
    result[2] = 6*x0*x0;
    result[1] = -4*x0*x0*x0;
    result[0] = x0*x0*x0*x0;
}

static inline void
dual2_eivals(dual2 a, complex result[4])
{
    double x0 = a.data[0];
    result[0] = complex_make(x0, 0);
    result[1] = complex_make(x0, 0);
    result[2] = complex_make(x0, 0);
    result[3] = complex_make(x0, 0);
}

#endif /* ALGPY_DUAL2_H */