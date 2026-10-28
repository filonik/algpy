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

#include "../../multicomplex/complex1.h"
#include "../../multicomplex/splitcomplex1.h"
#include "../../multicomplex/dual1.h"

#include "../../../fields/complex.h"
#include "../../../fields/double.h"
#include "../../../vectorspaces/vec1.h"
#include "../../../vectorspaces/vec2.h"

typedef vec2 multicomplex1;

static inline multicomplex1
multicomplex1_make(double x, double x0)
{
    multicomplex1 r = {{x, x0}};
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

/* Normalize e to square +1 or -1; zero uses the real derivative.
 * cf is the complex continuation of f on its principal branch. */
static inline multicomplex1
multicomplex1_apply(vec1 k, multicomplex1 a, double (*f)(double), double (*deriv1)(double), complex (*cf)(complex))
{
    if (k.data[0] == 0) {
        return dual1_apply(a, f, deriv1);
    }
    double s = sqrt(fabs(k.data[0]));
    complex value = complex_make(a.data[0], s*a.data[1]);
    if (k.data[0] > 0) {
        value = splitcomplex1_apply(value, f);
    } else if (k.data[0] < 0) {
        value = cf(value);
    } else {
        return multicomplex1_make(NAN, NAN);
    }
    return multicomplex1_make(value.data[0], value.data[1]/s);
}

static inline multicomplex1
multicomplex1_exp(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_exp(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_exp(a);
    }
    if (k.data[0] == 0) {
        return dual1_exp(a);
    }
    return multicomplex1_apply(k, a, double_exp, double_exp_deriv1, complex_exp);
}

static inline multicomplex1
multicomplex1_log(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_log(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_log(a);
    }
    if (k.data[0] == 0) {
        return dual1_log(a);
    }
    return multicomplex1_apply(k, a, double_log, double_log_deriv1, complex_log);
}

static inline multicomplex1
multicomplex1_pow(vec1 k, multicomplex1 a, multicomplex1 b)
{
    if (k.data[0] == -1) {
        return complex1_pow(a, b);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_pow(a, b);
    }
    if (k.data[0] == 0) {
        return dual1_pow(a, b);
    }
    double s = sqrt(fabs(k.data[0]));
    complex x = complex_make(a.data[0], s*a.data[1]);
    complex y = complex_make(b.data[0], s*b.data[1]);
    complex result;
    if (k.data[0] > 0) {
        result = splitcomplex1_pow(x, y);
    } else if (k.data[0] < 0) {
        result = complex_pow(x, y);
    } else {
        return multicomplex1_make(NAN, NAN);
    }
    return multicomplex1_make(result.data[0], result.data[1]/s);
}

static inline multicomplex1
multicomplex1_radd(double a, multicomplex1 b) {
  return vec2_add(multicomplex1_from_re(a), b);
}

static inline multicomplex1
multicomplex1_ladd(multicomplex1 a, double b) {
  return vec2_add(a, multicomplex1_from_re(b));
}

static inline multicomplex1
multicomplex1_rsub(double a, multicomplex1 b) {
  return vec2_sub(multicomplex1_from_re(a), b);
}

static inline multicomplex1
multicomplex1_lsub(multicomplex1 a, double b) {
  return vec2_sub(a, multicomplex1_from_re(b));
}

static inline multicomplex1
multicomplex1_rmul(vec1 k, double a, multicomplex1 b) {
  return multicomplex1_mul(k, multicomplex1_from_re(a), b);
}

static inline multicomplex1
multicomplex1_lmul(vec1 k, multicomplex1 a, double b) {
  return multicomplex1_mul(k, a, multicomplex1_from_re(b));
}

static inline multicomplex1
multicomplex1_rdiv(vec1 k, double a, multicomplex1 b) {
  return multicomplex1_div(k, multicomplex1_from_re(a), b);
}

static inline multicomplex1
multicomplex1_ldiv(vec1 k, multicomplex1 a, double b) {
  return multicomplex1_div(k, a, multicomplex1_from_re(b));
}

static inline multicomplex1
multicomplex1_rpow(vec1 k, double a, multicomplex1 b) {
  return multicomplex1_pow(k, multicomplex1_from_re(a), b);
}

static inline multicomplex1
multicomplex1_lpow(vec1 k, multicomplex1 a, double b) {
  return multicomplex1_pow(k, a, multicomplex1_from_re(b));
}

static inline multicomplex1
multicomplex1_mix(vec1 k, multicomplex1 a, multicomplex1 b, double t)
{
    multicomplex1 delta = multicomplex1_mul(k, multicomplex1_inv(k, a), b);
    return multicomplex1_mul(k, a, multicomplex1_lpow(k, delta, t));
}

static inline multicomplex1
multicomplex1_expm1(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_expm1(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_expm1(a);
    }
    if (k.data[0] == 0) {
        return dual1_expm1(a);
    }
    return multicomplex1_apply(k, a, double_expm1, double_expm1_deriv1, complex_expm1);
}

static inline multicomplex1
multicomplex1_log1p(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_log1p(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_log1p(a);
    }
    if (k.data[0] == 0) {
        return dual1_log1p(a);
    }
    return multicomplex1_apply(k, a, double_log1p, double_log1p_deriv1, complex_log1p);
}

static inline multicomplex1
multicomplex1_sqrt(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_sqrt(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_sqrt(a);
    }
    if (k.data[0] == 0) {
        return dual1_sqrt(a);
    }
    return multicomplex1_apply(k, a, double_sqrt, double_sqrt_deriv1, complex_sqrt);
}

static inline multicomplex1
multicomplex1_cos(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_cos(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_cos(a);
    }
    if (k.data[0] == 0) {
        return dual1_cos(a);
    }
    return multicomplex1_apply(k, a, double_cos, double_cos_deriv1, complex_cos);
}

static inline multicomplex1
multicomplex1_sin(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_sin(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_sin(a);
    }
    if (k.data[0] == 0) {
        return dual1_sin(a);
    }
    return multicomplex1_apply(k, a, double_sin, double_sin_deriv1, complex_sin);
}

static inline multicomplex1
multicomplex1_tan(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_tan(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_tan(a);
    }
    if (k.data[0] == 0) {
        return dual1_tan(a);
    }
    return multicomplex1_apply(k, a, double_tan, double_tan_deriv1, complex_tan);
}

static inline multicomplex1
multicomplex1_acos(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_acos(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_acos(a);
    }
    if (k.data[0] == 0) {
        return dual1_acos(a);
    }
    return multicomplex1_apply(k, a, double_acos, double_acos_deriv1, complex_acos);
}

static inline multicomplex1
multicomplex1_asin(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_asin(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_asin(a);
    }
    if (k.data[0] == 0) {
        return dual1_asin(a);
    }
    return multicomplex1_apply(k, a, double_asin, double_asin_deriv1, complex_asin);
}

static inline multicomplex1
multicomplex1_atan(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_atan(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_atan(a);
    }
    if (k.data[0] == 0) {
        return dual1_atan(a);
    }
    return multicomplex1_apply(k, a, double_atan, double_atan_deriv1, complex_atan);
}

static inline multicomplex1
multicomplex1_cosh(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_cosh(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_cosh(a);
    }
    if (k.data[0] == 0) {
        return dual1_cosh(a);
    }
    return multicomplex1_apply(k, a, double_cosh, double_cosh_deriv1, complex_cosh);
}

static inline multicomplex1
multicomplex1_sinh(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_sinh(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_sinh(a);
    }
    if (k.data[0] == 0) {
        return dual1_sinh(a);
    }
    return multicomplex1_apply(k, a, double_sinh, double_sinh_deriv1, complex_sinh);
}

static inline multicomplex1
multicomplex1_tanh(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_tanh(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_tanh(a);
    }
    if (k.data[0] == 0) {
        return dual1_tanh(a);
    }
    return multicomplex1_apply(k, a, double_tanh, double_tanh_deriv1, complex_tanh);
}

static inline multicomplex1
multicomplex1_acosh(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_acosh(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_acosh(a);
    }
    if (k.data[0] == 0) {
        return dual1_acosh(a);
    }
    return multicomplex1_apply(k, a, double_acosh, double_acosh_deriv1, complex_acosh);
}

static inline multicomplex1
multicomplex1_asinh(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_asinh(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_asinh(a);
    }
    if (k.data[0] == 0) {
        return dual1_asinh(a);
    }
    return multicomplex1_apply(k, a, double_asinh, double_asinh_deriv1, complex_asinh);
}

static inline multicomplex1
multicomplex1_atanh(vec1 k, multicomplex1 a)
{
    if (k.data[0] == -1) {
        return complex1_atanh(a);
    }
    if (k.data[0] == 1) {
        return splitcomplex1_atanh(a);
    }
    if (k.data[0] == 0) {
        return dual1_atanh(a);
    }
    return multicomplex1_apply(k, a, double_atanh, double_atanh_deriv1, complex_atanh);
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
multicomplex1_eivals(vec1 k, multicomplex1 a, complex result[2])
{
    complex s = complex_sqrt_real(k.data[0]);
    complex t = complex_make(a.data[1], 0);
    complex d = complex_mul(t, s);
    result[0] = complex_make(a.data[0] - d.data[0], -d.data[1]);
    result[1] = complex_make(a.data[0] + d.data[0], +d.data[1]);
}

/* Common algebra API. */
static inline multicomplex1
multicomplex1_zero(void)
{
    return vec2_zero();
}

static inline int
multicomplex1_is_real(multicomplex1 a)
{
    return vec2_is_real(a);
}

static inline multicomplex1
multicomplex1_pos(multicomplex1 a)
{
    return vec2_pos(a);
}

static inline multicomplex1
multicomplex1_neg(multicomplex1 a)
{
    return vec2_neg(a);
}

static inline multicomplex1
multicomplex1_add(multicomplex1 a, multicomplex1 b)
{
    return vec2_add(a, b);
}

static inline multicomplex1
multicomplex1_sub(multicomplex1 a, multicomplex1 b)
{
    return vec2_sub(a, b);
}

static inline int
multicomplex1_equal(multicomplex1 a, multicomplex1 b)
{
    return vec2_all_eq(a, b);
}

static inline int
multicomplex1_not_equal(multicomplex1 a, multicomplex1 b)
{
    return !multicomplex1_equal(a, b);
}

static inline multicomplex1
multicomplex1_square(vec1 k, multicomplex1 a)
{
    return multicomplex1_mul(k, a, a);
}

static inline complex
multicomplex1_norm(vec1 k, multicomplex1 a)
{
    return complex_root_real(multicomplex1_norm2(k, a), 2);
}

/* Derivatives through order 3, for nilpotent extensions at higher levels. */
static inline multicomplex1
multicomplex1_exp_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_exp_deriv1, double_exp_deriv2, complex_exp_deriv1);
}

static inline multicomplex1
multicomplex1_exp_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_exp_deriv2, double_exp_deriv3, complex_exp_deriv2);
}

static inline multicomplex1
multicomplex1_exp_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_exp_deriv3, double_exp_deriv4, complex_exp_deriv3);
}

static inline multicomplex1
multicomplex1_expm1_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_expm1_deriv1, double_expm1_deriv2, complex_expm1_deriv1);
}

static inline multicomplex1
multicomplex1_expm1_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_expm1_deriv2, double_expm1_deriv3, complex_expm1_deriv2);
}

static inline multicomplex1
multicomplex1_expm1_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_expm1_deriv3, double_expm1_deriv4, complex_expm1_deriv3);
}

static inline multicomplex1
multicomplex1_log_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_log_deriv1, double_log_deriv2, complex_log_deriv1);
}

static inline multicomplex1
multicomplex1_log_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_log_deriv2, double_log_deriv3, complex_log_deriv2);
}

static inline multicomplex1
multicomplex1_log_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_log_deriv3, double_log_deriv4, complex_log_deriv3);
}

static inline multicomplex1
multicomplex1_log1p_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_log1p_deriv1, double_log1p_deriv2, complex_log1p_deriv1);
}

static inline multicomplex1
multicomplex1_log1p_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_log1p_deriv2, double_log1p_deriv3, complex_log1p_deriv2);
}

static inline multicomplex1
multicomplex1_log1p_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_log1p_deriv3, double_log1p_deriv4, complex_log1p_deriv3);
}

static inline multicomplex1
multicomplex1_inv_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_inv_deriv1, double_inv_deriv2, complex_inv_deriv1);
}

static inline multicomplex1
multicomplex1_inv_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_inv_deriv2, double_inv_deriv3, complex_inv_deriv2);
}

static inline multicomplex1
multicomplex1_inv_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_inv_deriv3, double_inv_deriv4, complex_inv_deriv3);
}

static inline multicomplex1
multicomplex1_square_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_square_deriv1, double_square_deriv2, complex_square_deriv1);
}

static inline multicomplex1
multicomplex1_square_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_square_deriv2, double_square_deriv3, complex_square_deriv2);
}

static inline multicomplex1
multicomplex1_square_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_square_deriv3, double_square_deriv4, complex_square_deriv3);
}

static inline multicomplex1
multicomplex1_sqrt_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sqrt_deriv1, double_sqrt_deriv2, complex_sqrt_deriv1);
}

static inline multicomplex1
multicomplex1_sqrt_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt_deriv2);
}

static inline multicomplex1
multicomplex1_sqrt_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sqrt_deriv3, double_sqrt_deriv4, complex_sqrt_deriv3);
}

static inline multicomplex1
multicomplex1_sin_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sin_deriv1, double_sin_deriv2, complex_sin_deriv1);
}

static inline multicomplex1
multicomplex1_sin_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sin_deriv2, double_sin_deriv3, complex_sin_deriv2);
}

static inline multicomplex1
multicomplex1_sin_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sin_deriv3, double_sin_deriv4, complex_sin_deriv3);
}

static inline multicomplex1
multicomplex1_cos_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_cos_deriv1, double_cos_deriv2, complex_cos_deriv1);
}

static inline multicomplex1
multicomplex1_cos_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_cos_deriv2, double_cos_deriv3, complex_cos_deriv2);
}

static inline multicomplex1
multicomplex1_cos_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_cos_deriv3, double_cos_deriv4, complex_cos_deriv3);
}

static inline multicomplex1
multicomplex1_tan_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_tan_deriv1, double_tan_deriv2, complex_tan_deriv1);
}

static inline multicomplex1
multicomplex1_tan_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_tan_deriv2, double_tan_deriv3, complex_tan_deriv2);
}

static inline multicomplex1
multicomplex1_tan_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_tan_deriv3, double_tan_deriv4, complex_tan_deriv3);
}

static inline multicomplex1
multicomplex1_asin_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_asin_deriv1, double_asin_deriv2, complex_asin_deriv1);
}

static inline multicomplex1
multicomplex1_asin_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_asin_deriv2, double_asin_deriv3, complex_asin_deriv2);
}

static inline multicomplex1
multicomplex1_asin_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_asin_deriv3, double_asin_deriv4, complex_asin_deriv3);
}

static inline multicomplex1
multicomplex1_acos_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_acos_deriv1, double_acos_deriv2, complex_acos_deriv1);
}

static inline multicomplex1
multicomplex1_acos_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_acos_deriv2, double_acos_deriv3, complex_acos_deriv2);
}

static inline multicomplex1
multicomplex1_acos_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_acos_deriv3, double_acos_deriv4, complex_acos_deriv3);
}

static inline multicomplex1
multicomplex1_atan_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_atan_deriv1, double_atan_deriv2, complex_atan_deriv1);
}

static inline multicomplex1
multicomplex1_atan_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_atan_deriv2, double_atan_deriv3, complex_atan_deriv2);
}

static inline multicomplex1
multicomplex1_atan_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_atan_deriv3, double_atan_deriv4, complex_atan_deriv3);
}

static inline multicomplex1
multicomplex1_sinh_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sinh_deriv1, double_sinh_deriv2, complex_sinh_deriv1);
}

static inline multicomplex1
multicomplex1_sinh_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sinh_deriv2, double_sinh_deriv3, complex_sinh_deriv2);
}

static inline multicomplex1
multicomplex1_sinh_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_sinh_deriv3, double_sinh_deriv4, complex_sinh_deriv3);
}

static inline multicomplex1
multicomplex1_cosh_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_cosh_deriv1, double_cosh_deriv2, complex_cosh_deriv1);
}

static inline multicomplex1
multicomplex1_cosh_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_cosh_deriv2, double_cosh_deriv3, complex_cosh_deriv2);
}

static inline multicomplex1
multicomplex1_cosh_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_cosh_deriv3, double_cosh_deriv4, complex_cosh_deriv3);
}

static inline multicomplex1
multicomplex1_tanh_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_tanh_deriv1, double_tanh_deriv2, complex_tanh_deriv1);
}

static inline multicomplex1
multicomplex1_tanh_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_tanh_deriv2, double_tanh_deriv3, complex_tanh_deriv2);
}

static inline multicomplex1
multicomplex1_tanh_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_tanh_deriv3, double_tanh_deriv4, complex_tanh_deriv3);
}

static inline multicomplex1
multicomplex1_asinh_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_asinh_deriv1, double_asinh_deriv2, complex_asinh_deriv1);
}

static inline multicomplex1
multicomplex1_asinh_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_asinh_deriv2, double_asinh_deriv3, complex_asinh_deriv2);
}

static inline multicomplex1
multicomplex1_asinh_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_asinh_deriv3, double_asinh_deriv4, complex_asinh_deriv3);
}

static inline multicomplex1
multicomplex1_acosh_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_acosh_deriv1, double_acosh_deriv2, complex_acosh_deriv1);
}

static inline multicomplex1
multicomplex1_acosh_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_acosh_deriv2, double_acosh_deriv3, complex_acosh_deriv2);
}

static inline multicomplex1
multicomplex1_acosh_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_acosh_deriv3, double_acosh_deriv4, complex_acosh_deriv3);
}

static inline multicomplex1
multicomplex1_atanh_deriv1(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_atanh_deriv1, double_atanh_deriv2, complex_atanh_deriv1);
}

static inline multicomplex1
multicomplex1_atanh_deriv2(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_atanh_deriv2, double_atanh_deriv3, complex_atanh_deriv2);
}

static inline multicomplex1
multicomplex1_atanh_deriv3(vec1 k, multicomplex1 a)
{
    return multicomplex1_apply(k, a, double_atanh_deriv3, double_atanh_deriv4, complex_atanh_deriv3);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
multicomplex1_lmat(vec1 k, multicomplex1 a, double result[ALGPY_MULTICOMPLEX1_N*ALGPY_MULTICOMPLEX1_N])
{
    multicomplex1 coeffs[ALGPY_MULTICOMPLEX1_N * ALGPY_MULTICOMPLEX1_N];
    int row, col, i;
    multicomplex1_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_MULTICOMPLEX1_N; row++) {
        for (col = 0; col < ALGPY_MULTICOMPLEX1_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_MULTICOMPLEX1_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_MULTICOMPLEX1_N + col].data[row];
            }
            result[row * ALGPY_MULTICOMPLEX1_N + col] = value;
        }
    }
}

static inline void
multicomplex1_rmat(vec1 k, multicomplex1 a, double result[ALGPY_MULTICOMPLEX1_N*ALGPY_MULTICOMPLEX1_N])
{
    multicomplex1 coeffs[ALGPY_MULTICOMPLEX1_N * ALGPY_MULTICOMPLEX1_N];
    int row, col, i;
    multicomplex1_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_MULTICOMPLEX1_N; row++) {
        for (col = 0; col < ALGPY_MULTICOMPLEX1_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_MULTICOMPLEX1_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_MULTICOMPLEX1_N + i].data[col];
            }
            result[row * ALGPY_MULTICOMPLEX1_N + col] = value;
        }
    }
}


#endif /* ALGPY_MULTICOMPLEX1_H */
