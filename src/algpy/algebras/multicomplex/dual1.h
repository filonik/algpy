/*
 * Arithmetic of dual numbers (x + y*e), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_DUAL1_H
#define ALGPY_DUAL1_H

#define ALGPY_DUAL1_N 2

#include <stdint.h>
#include "../../fields/double.h"
#include "../../fields/complex.h"
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
    double x = a.data[0];
    return dual1_make(1/x, (-a.data[1]/x)/x);
}

static inline dual1
dual1_div(dual1 a, dual1 b)
{
    return dual1_mul(a, dual1_inv(b));
}

/* f and deriv1 must be defined and differentiable at the real part.
 * Singular derivatives follow IEEE floating-point propagation. */
static inline dual1
dual1_apply(dual1 a, double (*f)(double), double (*deriv1)(double))
{
    return dual1_make(f(a.data[0]), deriv1(a.data[0])*a.data[1]);
}

static inline dual1
dual1_exp(dual1 a)
{
    return dual1_apply(a, double_exp, double_exp_deriv1);
}

static inline dual1
dual1_log(dual1 a)
{
    return dual1_apply(a, double_log, double_log_deriv1);
}

static inline dual1
dual1_sqrt(dual1 a)
{
    return dual1_apply(a, double_sqrt, double_sqrt_deriv1);
}

static inline dual1
dual1_cos(dual1 a)
{
    return dual1_apply(a, double_cos, double_cos_deriv1);
}

static inline dual1
dual1_sin(dual1 a)
{
    return dual1_apply(a, double_sin, double_sin_deriv1);
}

static inline dual1
dual1_tan(dual1 a)
{
    return dual1_apply(a, double_tan, double_tan_deriv1);
}

static inline dual1
dual1_acos(dual1 a)
{
    return dual1_apply(a, double_acos, double_acos_deriv1);
}

static inline dual1
dual1_asin(dual1 a)
{
    return dual1_apply(a, double_asin, double_asin_deriv1);
}

static inline dual1
dual1_atan(dual1 a)
{
    return dual1_apply(a, double_atan, double_atan_deriv1);
}

static inline dual1
dual1_cosh(dual1 a)
{
    return dual1_apply(a, double_cosh, double_cosh_deriv1);
}

static inline dual1
dual1_sinh(dual1 a)
{
    return dual1_apply(a, double_sinh, double_sinh_deriv1);
}

static inline dual1
dual1_tanh(dual1 a)
{
    return dual1_apply(a, double_tanh, double_tanh_deriv1);
}

static inline dual1
dual1_acosh(dual1 a)
{
    return dual1_apply(a, double_acosh, double_acosh_deriv1);
}

static inline dual1
dual1_asinh(dual1 a)
{
    return dual1_apply(a, double_asinh, double_asinh_deriv1);
}

static inline dual1
dual1_atanh(dual1 a)
{
    return dual1_apply(a, double_atanh, double_atanh_deriv1);
}

/* Real analytic power: requires a positive real part for log(a). */
/* Integer powers by squaring; negative exponents invert the base once. */
static inline dual1
dual1_ipow(dual1 a, int64_t b)
{
    dual1 result = dual1_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return dual1_mul(a, a);
    }
    if (b == -1) {
        return dual1_inv(a);
    }
    if (b < 0) {
        a = dual1_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = dual1_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = dual1_mul(a, a);
        }
    }
    return result;
}

static inline dual1
dual1_rpow(dual1 a, dual1 b)
{
    return dual1_exp(dual1_mul(dual1_log(a), b));
}

static inline dual1
dual1_pow(dual1 a, dual1 b)
{
    return dual1_rpow(a, b);
}

static inline dual1
dual1_lpow(dual1 a, dual1 b)
{
    return dual1_rpow(a, b);
}

static inline void
dual1_minpoly(dual1 a, double result[3])
{
    double x0 = a.data[0];
    result[2] = 1;
    result[1] = -2*x0;
    result[0] = x0*x0;
}

static inline void
dual1_eivals(dual1 a, complex result[2])
{
    double x0 = a.data[0];
    result[0] = complex_make(x0, 0);
    result[1] = complex_make(x0, 0);
}

/* Common algebra API. */
static inline dual1
dual1_zero(void)
{
    return vec2_zero();
}

static inline int
dual1_is_real(dual1 a)
{
    return vec2_is_real(a);
}

static inline dual1
dual1_pos(dual1 a)
{
    return vec2_pos(a);
}

static inline dual1
dual1_neg(dual1 a)
{
    return vec2_neg(a);
}

static inline dual1
dual1_add(dual1 a, dual1 b)
{
    return vec2_add(a, b);
}

static inline dual1
dual1_sub(dual1 a, dual1 b)
{
    return vec2_sub(a, b);
}

static inline int
dual1_equal(dual1 a, dual1 b)
{
    return vec2_all_eq(a, b);
}

static inline int
dual1_not_equal(dual1 a, dual1 b)
{
    return !dual1_equal(a, b);
}

static inline dual1
dual1_square(dual1 a)
{
    return dual1_mul(a, a);
}

static inline dual1
dual1_radd(double a, dual1 b)
{
    return dual1_add(dual1_from_re(a), b);
}

static inline dual1
dual1_ladd(dual1 a, double b)
{
    return dual1_add(a, dual1_from_re(b));
}

static inline dual1
dual1_rsub(double a, dual1 b)
{
    return dual1_sub(dual1_from_re(a), b);
}

static inline dual1
dual1_lsub(dual1 a, double b)
{
    return dual1_sub(a, dual1_from_re(b));
}

static inline dual1
dual1_rmul(double a, dual1 b)
{
    return dual1_mul(dual1_from_re(a), b);
}

static inline dual1
dual1_lmul(dual1 a, double b)
{
    return dual1_mul(a, dual1_from_re(b));
}

static inline dual1
dual1_rdiv(double a, dual1 b)
{
    return dual1_div(dual1_from_re(a), b);
}

static inline dual1
dual1_ldiv(dual1 a, double b)
{
    return dual1_div(a, dual1_from_re(b));
}

static inline dual1
dual1_rpows(double a, dual1 b)
{
    return dual1_rpow(dual1_from_re(a), b);
}

static inline dual1
dual1_lpows(dual1 a, double b)
{
    return dual1_rpow(a, dual1_from_re(b));
}

static inline dual1
dual1_expm1(dual1 a)
{
    return dual1_apply(a, double_expm1, double_expm1_deriv1);
}

static inline dual1
dual1_log1p(dual1 a)
{
    return dual1_apply(a, double_log1p, double_log1p_deriv1);
}

static inline dual1
dual1_mix(dual1 a, dual1 b, double t)
{
    dual1 delta = dual1_mul(dual1_inv(a), b);
    return dual1_mul(a, dual1_lpows(delta, t));
}

static inline double
dual1_norm(dual1 a)
{
    return fabs(a.data[0]);
}

/* Derivatives through order 3; each dual level uses one more field derivative. */
static inline dual1
dual1_exp_deriv1(dual1 a)
{
    return dual1_apply(a, double_exp_deriv1, double_exp_deriv2);
}

static inline dual1
dual1_exp_deriv2(dual1 a)
{
    return dual1_apply(a, double_exp_deriv2, double_exp_deriv3);
}

static inline dual1
dual1_exp_deriv3(dual1 a)
{
    return dual1_apply(a, double_exp_deriv3, double_exp_deriv4);
}

static inline dual1
dual1_expm1_deriv1(dual1 a)
{
    return dual1_apply(a, double_expm1_deriv1, double_expm1_deriv2);
}

static inline dual1
dual1_expm1_deriv2(dual1 a)
{
    return dual1_apply(a, double_expm1_deriv2, double_expm1_deriv3);
}

static inline dual1
dual1_expm1_deriv3(dual1 a)
{
    return dual1_apply(a, double_expm1_deriv3, double_expm1_deriv4);
}

static inline dual1
dual1_log_deriv1(dual1 a)
{
    return dual1_apply(a, double_log_deriv1, double_log_deriv2);
}

static inline dual1
dual1_log_deriv2(dual1 a)
{
    return dual1_apply(a, double_log_deriv2, double_log_deriv3);
}

static inline dual1
dual1_log_deriv3(dual1 a)
{
    return dual1_apply(a, double_log_deriv3, double_log_deriv4);
}

static inline dual1
dual1_log1p_deriv1(dual1 a)
{
    return dual1_apply(a, double_log1p_deriv1, double_log1p_deriv2);
}

static inline dual1
dual1_log1p_deriv2(dual1 a)
{
    return dual1_apply(a, double_log1p_deriv2, double_log1p_deriv3);
}

static inline dual1
dual1_log1p_deriv3(dual1 a)
{
    return dual1_apply(a, double_log1p_deriv3, double_log1p_deriv4);
}

static inline dual1
dual1_inv_deriv1(dual1 a)
{
    return dual1_apply(a, double_inv_deriv1, double_inv_deriv2);
}

static inline dual1
dual1_inv_deriv2(dual1 a)
{
    return dual1_apply(a, double_inv_deriv2, double_inv_deriv3);
}

static inline dual1
dual1_inv_deriv3(dual1 a)
{
    return dual1_apply(a, double_inv_deriv3, double_inv_deriv4);
}

static inline dual1
dual1_square_deriv1(dual1 a)
{
    return dual1_apply(a, double_square_deriv1, double_square_deriv2);
}

static inline dual1
dual1_square_deriv2(dual1 a)
{
    return dual1_apply(a, double_square_deriv2, double_square_deriv3);
}

static inline dual1
dual1_square_deriv3(dual1 a)
{
    return dual1_apply(a, double_square_deriv3, double_square_deriv4);
}

static inline dual1
dual1_sqrt_deriv1(dual1 a)
{
    return dual1_apply(a, double_sqrt_deriv1, double_sqrt_deriv2);
}

static inline dual1
dual1_sqrt_deriv2(dual1 a)
{
    return dual1_apply(a, double_sqrt_deriv2, double_sqrt_deriv3);
}

static inline dual1
dual1_sqrt_deriv3(dual1 a)
{
    return dual1_apply(a, double_sqrt_deriv3, double_sqrt_deriv4);
}

static inline dual1
dual1_sin_deriv1(dual1 a)
{
    return dual1_apply(a, double_sin_deriv1, double_sin_deriv2);
}

static inline dual1
dual1_sin_deriv2(dual1 a)
{
    return dual1_apply(a, double_sin_deriv2, double_sin_deriv3);
}

static inline dual1
dual1_sin_deriv3(dual1 a)
{
    return dual1_apply(a, double_sin_deriv3, double_sin_deriv4);
}

static inline dual1
dual1_cos_deriv1(dual1 a)
{
    return dual1_apply(a, double_cos_deriv1, double_cos_deriv2);
}

static inline dual1
dual1_cos_deriv2(dual1 a)
{
    return dual1_apply(a, double_cos_deriv2, double_cos_deriv3);
}

static inline dual1
dual1_cos_deriv3(dual1 a)
{
    return dual1_apply(a, double_cos_deriv3, double_cos_deriv4);
}

static inline dual1
dual1_tan_deriv1(dual1 a)
{
    return dual1_apply(a, double_tan_deriv1, double_tan_deriv2);
}

static inline dual1
dual1_tan_deriv2(dual1 a)
{
    return dual1_apply(a, double_tan_deriv2, double_tan_deriv3);
}

static inline dual1
dual1_tan_deriv3(dual1 a)
{
    return dual1_apply(a, double_tan_deriv3, double_tan_deriv4);
}

static inline dual1
dual1_asin_deriv1(dual1 a)
{
    return dual1_apply(a, double_asin_deriv1, double_asin_deriv2);
}

static inline dual1
dual1_asin_deriv2(dual1 a)
{
    return dual1_apply(a, double_asin_deriv2, double_asin_deriv3);
}

static inline dual1
dual1_asin_deriv3(dual1 a)
{
    return dual1_apply(a, double_asin_deriv3, double_asin_deriv4);
}

static inline dual1
dual1_acos_deriv1(dual1 a)
{
    return dual1_apply(a, double_acos_deriv1, double_acos_deriv2);
}

static inline dual1
dual1_acos_deriv2(dual1 a)
{
    return dual1_apply(a, double_acos_deriv2, double_acos_deriv3);
}

static inline dual1
dual1_acos_deriv3(dual1 a)
{
    return dual1_apply(a, double_acos_deriv3, double_acos_deriv4);
}

static inline dual1
dual1_atan_deriv1(dual1 a)
{
    return dual1_apply(a, double_atan_deriv1, double_atan_deriv2);
}

static inline dual1
dual1_atan_deriv2(dual1 a)
{
    return dual1_apply(a, double_atan_deriv2, double_atan_deriv3);
}

static inline dual1
dual1_atan_deriv3(dual1 a)
{
    return dual1_apply(a, double_atan_deriv3, double_atan_deriv4);
}

static inline dual1
dual1_sinh_deriv1(dual1 a)
{
    return dual1_apply(a, double_sinh_deriv1, double_sinh_deriv2);
}

static inline dual1
dual1_sinh_deriv2(dual1 a)
{
    return dual1_apply(a, double_sinh_deriv2, double_sinh_deriv3);
}

static inline dual1
dual1_sinh_deriv3(dual1 a)
{
    return dual1_apply(a, double_sinh_deriv3, double_sinh_deriv4);
}

static inline dual1
dual1_cosh_deriv1(dual1 a)
{
    return dual1_apply(a, double_cosh_deriv1, double_cosh_deriv2);
}

static inline dual1
dual1_cosh_deriv2(dual1 a)
{
    return dual1_apply(a, double_cosh_deriv2, double_cosh_deriv3);
}

static inline dual1
dual1_cosh_deriv3(dual1 a)
{
    return dual1_apply(a, double_cosh_deriv3, double_cosh_deriv4);
}

static inline dual1
dual1_tanh_deriv1(dual1 a)
{
    return dual1_apply(a, double_tanh_deriv1, double_tanh_deriv2);
}

static inline dual1
dual1_tanh_deriv2(dual1 a)
{
    return dual1_apply(a, double_tanh_deriv2, double_tanh_deriv3);
}

static inline dual1
dual1_tanh_deriv3(dual1 a)
{
    return dual1_apply(a, double_tanh_deriv3, double_tanh_deriv4);
}

static inline dual1
dual1_asinh_deriv1(dual1 a)
{
    return dual1_apply(a, double_asinh_deriv1, double_asinh_deriv2);
}

static inline dual1
dual1_asinh_deriv2(dual1 a)
{
    return dual1_apply(a, double_asinh_deriv2, double_asinh_deriv3);
}

static inline dual1
dual1_asinh_deriv3(dual1 a)
{
    return dual1_apply(a, double_asinh_deriv3, double_asinh_deriv4);
}

static inline dual1
dual1_acosh_deriv1(dual1 a)
{
    return dual1_apply(a, double_acosh_deriv1, double_acosh_deriv2);
}

static inline dual1
dual1_acosh_deriv2(dual1 a)
{
    return dual1_apply(a, double_acosh_deriv2, double_acosh_deriv3);
}

static inline dual1
dual1_acosh_deriv3(dual1 a)
{
    return dual1_apply(a, double_acosh_deriv3, double_acosh_deriv4);
}

static inline dual1
dual1_atanh_deriv1(dual1 a)
{
    return dual1_apply(a, double_atanh_deriv1, double_atanh_deriv2);
}

static inline dual1
dual1_atanh_deriv2(dual1 a)
{
    return dual1_apply(a, double_atanh_deriv2, double_atanh_deriv3);
}

static inline dual1
dual1_atanh_deriv3(dual1 a)
{
    return dual1_apply(a, double_atanh_deriv3, double_atanh_deriv4);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
dual1_lmat(dual1 a, double result[ALGPY_DUAL1_N*ALGPY_DUAL1_N])
{
    dual1 basis = {{0}};
    dual1 product;
    int row, col;
    for (col = 0; col < ALGPY_DUAL1_N; col++) {
        basis.data[col] = 1;
        product = dual1_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_DUAL1_N; row++) {
            result[row * ALGPY_DUAL1_N + col] = product.data[row];
        }
    }
}

static inline void
dual1_rmat(dual1 a, double result[ALGPY_DUAL1_N*ALGPY_DUAL1_N])
{
    dual1 basis = {{0}};
    dual1 product;
    int row, col;
    for (row = 0; row < ALGPY_DUAL1_N; row++) {
        basis.data[row] = 1;
        product = dual1_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_DUAL1_N; col++) {
            result[row * ALGPY_DUAL1_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_DUAL1_H */
