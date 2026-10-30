/*
 * Arithmetic of bi-dual numbers x + x0*e0 + x1*e1 + x01*e0*e1 with
 * e0*e0 = e1*e1 = 0, stored as a vec4.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_DUAL2_H
#define ALGPY_DUAL2_H

#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../vectorspaces/vec4.h"
#include "dual1.h"

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

static inline double
dual2_norm(dual2 a)
{
    return fabs(a.data[0]);
}

static inline dual2
dual2_inv(dual2 a)
{
    dual1 u = dual1_make(a.data[0], a.data[1]);
    dual1 v = dual1_make(a.data[2], a.data[3]);
    dual1 ru = dual1_inv(u);
    dual1 rv = vec2_neg(dual1_mul(dual1_mul(v, ru), ru));
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

/* With a = u + v*e1 and e1^2 = 0, f(a) = f(u) + v*deriv1(u)*e1.
 * deriv1 itself is dual-valued, so it carries the mixed second derivative.
 * Callbacks must be defined and differentiable at the real part.
 */
static inline dual2
dual2_apply(dual2 a, dual1 (*f)(dual1), dual1 (*deriv1)(dual1))
{
    dual1 u = dual1_make(a.data[0], a.data[1]);
    dual1 v = dual1_make(a.data[2], a.data[3]);
    dual1 value = f(u);
    dual1 tangent = dual1_mul(v, deriv1(u));
    return dual2_make(value.data[0], value.data[1], tangent.data[0], tangent.data[1]);
}

static inline dual2
dual2_exp(dual2 a)
{
    return dual2_apply(a, dual1_exp, dual1_exp_deriv1);
}

static inline dual2
dual2_log(dual2 a)
{
    return dual2_apply(a, dual1_log, dual1_log_deriv1);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline dual2
dual2_ipow(dual2 a, int64_t b)
{
    dual2 result = dual2_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return dual2_mul(a, a);
    }
    if (b == -1) {
        return dual2_inv(a);
    }
    if (b < 0) {
        a = dual2_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = dual2_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = dual2_mul(a, a);
        }
    }
    return result;
}

static inline dual2
dual2_rpow(dual2 a, dual2 b)
{
    return dual2_exp(dual2_mul(dual2_log(a), b));
}

static inline dual2
dual2_pow(dual2 a, dual2 b)
{
    return dual2_rpow(a, b);
}

static inline dual2
dual2_lpow(dual2 a, dual2 b)
{
    return dual2_rpow(a, b);
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
dual2_rpows(double a, dual2 b) {
  return dual2_rpow(dual2_from_re(a), b);
}

static inline dual2
dual2_lpows(dual2 a, double b) {
  return dual2_rpow(a, dual2_from_re(b));
}

static inline dual2
dual2_mix(dual2 a, dual2 b, double t)
{
    dual2 delta = dual2_mul(dual2_inv(a), b);
    return dual2_mul(a, dual2_lpows(delta, t));
}

static inline dual2
dual2_cos(dual2 a)
{
    return dual2_apply(a, dual1_cos, dual1_cos_deriv1);
}

static inline dual2
dual2_sin(dual2 a)
{
    return dual2_apply(a, dual1_sin, dual1_sin_deriv1);
}

static inline dual2
dual2_tan(dual2 a)
{
    return dual2_apply(a, dual1_tan, dual1_tan_deriv1);
}

static inline dual2
dual2_cosh(dual2 a)
{
    return dual2_apply(a, dual1_cosh, dual1_cosh_deriv1);
}

static inline dual2
dual2_sinh(dual2 a)
{
    return dual2_apply(a, dual1_sinh, dual1_sinh_deriv1);
}

static inline dual2
dual2_tanh(dual2 a)
{
    return dual2_apply(a, dual1_tanh, dual1_tanh_deriv1);
}

static inline dual2
dual2_sqrt(dual2 a)
{
    return dual2_apply(a, dual1_sqrt, dual1_sqrt_deriv1);
}

static inline dual2
dual2_acos(dual2 a)
{
    return dual2_apply(a, dual1_acos, dual1_acos_deriv1);
}

static inline dual2
dual2_asin(dual2 a)
{
    return dual2_apply(a, dual1_asin, dual1_asin_deriv1);
}

static inline dual2
dual2_atan(dual2 a)
{
    return dual2_apply(a, dual1_atan, dual1_atan_deriv1);
}

static inline dual2
dual2_acosh(dual2 a)
{
    return dual2_apply(a, dual1_acosh, dual1_acosh_deriv1);
}

static inline dual2
dual2_asinh(dual2 a)
{
    return dual2_apply(a, dual1_asinh, dual1_asinh_deriv1);
}

static inline dual2
dual2_atanh(dual2 a)
{
    return dual2_apply(a, dual1_atanh, dual1_atanh_deriv1);
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

/* Common algebra API. */
static inline dual2
dual2_zero(void)
{
    return vec4_zero();
}

static inline int
dual2_is_real(dual2 a)
{
    return vec4_is_real(a);
}

static inline dual2
dual2_pos(dual2 a)
{
    return vec4_pos(a);
}

static inline dual2
dual2_neg(dual2 a)
{
    return vec4_neg(a);
}

static inline dual2
dual2_add(dual2 a, dual2 b)
{
    return vec4_add(a, b);
}

static inline dual2
dual2_sub(dual2 a, dual2 b)
{
    return vec4_sub(a, b);
}

static inline int
dual2_equal(dual2 a, dual2 b)
{
    return vec4_all_eq(a, b);
}

static inline int
dual2_not_equal(dual2 a, dual2 b)
{
    return !dual2_equal(a, b);
}

static inline dual2
dual2_re(dual2 a)
{
    return dual2_from_re(a.data[0]);
}

static inline dual2
dual2_im(dual2 a)
{
    return dual2_sub(a, dual2_re(a));
}

static inline dual2
dual2_square(dual2 a)
{
    return dual2_mul(a, a);
}

static inline dual2
dual2_expm1(dual2 a)
{
    return dual2_apply(a, dual1_expm1, dual1_expm1_deriv1);
}

static inline dual2
dual2_log1p(dual2 a)
{
    return dual2_apply(a, dual1_log1p, dual1_log1p_deriv1);
}

/* Derivatives through order 2; each dual level uses one more field derivative. */
static inline dual2
dual2_exp_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_exp_deriv1, dual1_exp_deriv2);
}

static inline dual2
dual2_exp_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_exp_deriv2, dual1_exp_deriv3);
}

static inline dual2
dual2_expm1_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_expm1_deriv1, dual1_expm1_deriv2);
}

static inline dual2
dual2_expm1_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_expm1_deriv2, dual1_expm1_deriv3);
}

static inline dual2
dual2_log_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_log_deriv1, dual1_log_deriv2);
}

static inline dual2
dual2_log_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_log_deriv2, dual1_log_deriv3);
}

static inline dual2
dual2_log1p_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_log1p_deriv1, dual1_log1p_deriv2);
}

static inline dual2
dual2_log1p_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_log1p_deriv2, dual1_log1p_deriv3);
}

static inline dual2
dual2_inv_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_inv_deriv1, dual1_inv_deriv2);
}

static inline dual2
dual2_inv_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_inv_deriv2, dual1_inv_deriv3);
}

static inline dual2
dual2_square_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_square_deriv1, dual1_square_deriv2);
}

static inline dual2
dual2_square_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_square_deriv2, dual1_square_deriv3);
}

static inline dual2
dual2_sqrt_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_sqrt_deriv1, dual1_sqrt_deriv2);
}

static inline dual2
dual2_sqrt_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_sqrt_deriv2, dual1_sqrt_deriv3);
}

static inline dual2
dual2_sin_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_sin_deriv1, dual1_sin_deriv2);
}

static inline dual2
dual2_sin_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_sin_deriv2, dual1_sin_deriv3);
}

static inline dual2
dual2_cos_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_cos_deriv1, dual1_cos_deriv2);
}

static inline dual2
dual2_cos_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_cos_deriv2, dual1_cos_deriv3);
}

static inline dual2
dual2_tan_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_tan_deriv1, dual1_tan_deriv2);
}

static inline dual2
dual2_tan_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_tan_deriv2, dual1_tan_deriv3);
}

static inline dual2
dual2_asin_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_asin_deriv1, dual1_asin_deriv2);
}

static inline dual2
dual2_asin_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_asin_deriv2, dual1_asin_deriv3);
}

static inline dual2
dual2_acos_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_acos_deriv1, dual1_acos_deriv2);
}

static inline dual2
dual2_acos_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_acos_deriv2, dual1_acos_deriv3);
}

static inline dual2
dual2_atan_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_atan_deriv1, dual1_atan_deriv2);
}

static inline dual2
dual2_atan_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_atan_deriv2, dual1_atan_deriv3);
}

static inline dual2
dual2_sinh_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_sinh_deriv1, dual1_sinh_deriv2);
}

static inline dual2
dual2_sinh_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_sinh_deriv2, dual1_sinh_deriv3);
}

static inline dual2
dual2_cosh_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_cosh_deriv1, dual1_cosh_deriv2);
}

static inline dual2
dual2_cosh_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_cosh_deriv2, dual1_cosh_deriv3);
}

static inline dual2
dual2_tanh_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_tanh_deriv1, dual1_tanh_deriv2);
}

static inline dual2
dual2_tanh_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_tanh_deriv2, dual1_tanh_deriv3);
}

static inline dual2
dual2_asinh_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_asinh_deriv1, dual1_asinh_deriv2);
}

static inline dual2
dual2_asinh_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_asinh_deriv2, dual1_asinh_deriv3);
}

static inline dual2
dual2_acosh_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_acosh_deriv1, dual1_acosh_deriv2);
}

static inline dual2
dual2_acosh_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_acosh_deriv2, dual1_acosh_deriv3);
}

static inline dual2
dual2_atanh_deriv1(dual2 a)
{
    return dual2_apply(a, dual1_atanh_deriv1, dual1_atanh_deriv2);
}

static inline dual2
dual2_atanh_deriv2(dual2 a)
{
    return dual2_apply(a, dual1_atanh_deriv2, dual1_atanh_deriv3);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
dual2_lmat(dual2 a, double result[ALGPY_DUAL2_N*ALGPY_DUAL2_N])
{
    dual2 basis = {{0}};
    dual2 product;
    int row, col;
    for (col = 0; col < ALGPY_DUAL2_N; col++) {
        basis.data[col] = 1;
        product = dual2_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_DUAL2_N; row++) {
            result[row * ALGPY_DUAL2_N + col] = product.data[row];
        }
    }
}

static inline void
dual2_rmat(dual2 a, double result[ALGPY_DUAL2_N*ALGPY_DUAL2_N])
{
    dual2 basis = {{0}};
    dual2 product;
    int row, col;
    for (row = 0; row < ALGPY_DUAL2_N; row++) {
        basis.data[row] = 1;
        product = dual2_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_DUAL2_N; col++) {
            result[row * ALGPY_DUAL2_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_DUAL2_H */
