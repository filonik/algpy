/*
 * Jets of order 5: truncated polynomials
 * x + x1*e + ... + x5*e^5 in R[e]/(e^6), stored as a vec6.
 * A function applied to a = x + e gives its Taylor coefficients, the
 * normalized derivatives f^(k)(x)/k!, k = 0..5.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_JET5_H
#define ALGPY_JET5_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../vectorspaces/vec6.h"

#define ALGPY_JET5_N 6

typedef vec6 jet5;

static inline jet5
jet5_make(const double coeffs[ALGPY_JET5_N])
{
    jet5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline jet5
jet5_from_re(double a)
{
    return jet5_make((const double[ALGPY_JET5_N]){a});
}

static inline double
jet5_to_re(jet5 a)
{
    return a.data[0];
}

static inline jet5
jet5_zero(void)
{
    return vec6_zero();
}

static inline jet5
jet5_one(void)
{
    return jet5_from_re(1);
}

static inline int
jet5_is_real(jet5 a)
{
    return vec6_is_real(a);
}

static inline jet5
jet5_re(jet5 a)
{
    return jet5_from_re(a.data[0]);
}

static inline jet5
jet5_im(jet5 a)
{
    jet5 result = a;
    result.data[0] = 0;
    return result;
}

static inline jet5
jet5_pos(jet5 a)
{
    return vec6_pos(a);
}

static inline jet5
jet5_neg(jet5 a)
{
    return vec6_neg(a);
}

static inline jet5
jet5_add(jet5 a, jet5 b)
{
    return vec6_add(a, b);
}

static inline jet5
jet5_sub(jet5 a, jet5 b)
{
    return vec6_sub(a, b);
}

static inline int
jet5_equal(jet5 a, jet5 b)
{
    return vec6_all_eq(a, b);
}

static inline int
jet5_not_equal(jet5 a, jet5 b)
{
    return !jet5_equal(a, b);
}

/* Truncated Cauchy product: coefficient k is the sum of x[i]*y[k-i]. */
static inline jet5
jet5_mul(jet5 a, jet5 b)
{
    const double *x = a.data, *y = b.data;
    return jet5_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] + x[1]*y[1] + x[2]*y[0],
        x[0]*y[3] + x[1]*y[2] + x[2]*y[1] + x[3]*y[0],
        x[0]*y[4] + x[1]*y[3] + x[2]*y[2] + x[3]*y[1] + x[4]*y[0],
        x[0]*y[5] + x[1]*y[4] + x[2]*y[3] + x[3]*y[2] + x[4]*y[1] + x[5]*y[0]});
}

static inline jet5
jet5_square(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){
        x[0]*x[0],
        2*x[0]*x[1],
        2*x[0]*x[2] + x[1]*x[1],
        2*x[0]*x[3] + 2*x[1]*x[2],
        2*x[0]*x[4] + 2*x[1]*x[3] + x[2]*x[2],
        2*x[0]*x[5] + 2*x[1]*x[4] + 2*x[2]*x[3]});
}

/* The automorphism e -> -e: coefficient k changes sign for odd k. */
static inline jet5
jet5_conjugate(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){x[0], -x[1], x[2], -x[3], x[4], -x[5]});
}

static inline void
jet5_basis(jet5 result[ALGPY_JET5_N])
{
    int i, j;
    for (i = 0; i < ALGPY_JET5_N; i++) {
        for (j = 0; j < ALGPY_JET5_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

/* The single generator e, at coefficient position 1. */
static inline void
jet5_gens(jet5 result[1])
{
    result[0] = jet5_make((const double[ALGPY_JET5_N]){[1] = 1});
}

static inline void
jet5_coeffs(jet5 result[ALGPY_JET5_N * ALGPY_JET5_N])
{
    jet5 basis[ALGPY_JET5_N];
    int i, j;
    jet5_basis(basis);
    for (i = 0; i < ALGPY_JET5_N; i++) {
        for (j = 0; j < ALGPY_JET5_N; j++) {
            result[i * ALGPY_JET5_N + j] = jet5_mul(basis[i], basis[j]);
        }
    }
}

/* Grade projections: grade k is the coefficient of e^k. */
static inline jet5
jet5_grade0(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){x[0], 0, 0, 0, 0, 0});
}

static inline jet5
jet5_grade1(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){0, x[1], 0, 0, 0, 0});
}

static inline jet5
jet5_grade2(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){0, 0, x[2], 0, 0, 0});
}

static inline jet5
jet5_grade3(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){0, 0, 0, x[3], 0, 0});
}

static inline jet5
jet5_grade4(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){0, 0, 0, 0, x[4], 0});
}

static inline jet5
jet5_grade5(jet5 a)
{
    const double *x = a.data;
    return jet5_make((const double[]){0, 0, 0, 0, 0, x[5]});
}

/* Reciprocal by the series recurrence y0 = 1/x0 and
 * y_k = -(x1*y_(k-1) + ... + xk*y0)/x0. Requires a nonzero real part. */
static inline jet5
jet5_inv(jet5 a)
{
    const double *x = a.data;
    jet5 result;
    double *y = result.data;
    int k, i;
    y[0] = 1/x[0];
    for (k = 1; k < ALGPY_JET5_N; k++) {
        double sum = 0;
        for (i = 1; i <= k; i++) {
            sum += x[i]*y[k - i];
        }
        y[k] = -sum/x[0];
    }
    return result;
}

static inline jet5
jet5_div(jet5 a, jet5 b)
{
    if (vec6_is_real(b)) {
        return vec6_div_scalar(a, b.data[0]);
    }
    return jet5_mul(a, jet5_inv(b));
}

static const double jet5_factorials[6] = {1., 1., 2., 6., 24., 120.};

/*
 * Taylor expansion about the real part x: with a = x + n, n nilpotent,
 * f(a) = sum over k of f^(k)(x)/k! * n^k, evaluated by Horner's rule.
 * derivs[k-1] is the k-th derivative; f must be 5 times
 * differentiable at x. Singular derivatives follow IEEE propagation.
 */
static inline jet5
jet5_apply(jet5 a, double (*f)(double), double (*const derivs[5])(double))
{
    double x = a.data[0];
    jet5 n = a, result = {{0}};
    int k;
    n.data[0] = 0;
    result.data[0] = derivs[4](x)/120.;
    for (k = 4; k >= 1; k--) {
        result = jet5_mul(result, n);
        result.data[0] += derivs[k - 1](x)/jet5_factorials[k];
    }
    result = jet5_mul(result, n);
    result.data[0] += f(x);
    return result;
}

static inline jet5
jet5_exp(jet5 a)
{
    static double (*const derivs[5])(double) = {double_exp_deriv1, double_exp_deriv2, double_exp_deriv3, double_exp_deriv4, double_exp_deriv5};
    return jet5_apply(a, double_exp, derivs);
}

static inline jet5
jet5_expm1(jet5 a)
{
    static double (*const derivs[5])(double) = {double_expm1_deriv1, double_expm1_deriv2, double_expm1_deriv3, double_expm1_deriv4, double_expm1_deriv5};
    return jet5_apply(a, double_expm1, derivs);
}

static inline jet5
jet5_log(jet5 a)
{
    static double (*const derivs[5])(double) = {double_log_deriv1, double_log_deriv2, double_log_deriv3, double_log_deriv4, double_log_deriv5};
    return jet5_apply(a, double_log, derivs);
}

static inline jet5
jet5_log1p(jet5 a)
{
    static double (*const derivs[5])(double) = {double_log1p_deriv1, double_log1p_deriv2, double_log1p_deriv3, double_log1p_deriv4, double_log1p_deriv5};
    return jet5_apply(a, double_log1p, derivs);
}

static inline jet5
jet5_sqrt(jet5 a)
{
    static double (*const derivs[5])(double) = {double_sqrt_deriv1, double_sqrt_deriv2, double_sqrt_deriv3, double_sqrt_deriv4, double_sqrt_deriv5};
    return jet5_apply(a, double_sqrt, derivs);
}

static inline jet5
jet5_sin(jet5 a)
{
    static double (*const derivs[5])(double) = {double_sin_deriv1, double_sin_deriv2, double_sin_deriv3, double_sin_deriv4, double_sin_deriv5};
    return jet5_apply(a, double_sin, derivs);
}

static inline jet5
jet5_cos(jet5 a)
{
    static double (*const derivs[5])(double) = {double_cos_deriv1, double_cos_deriv2, double_cos_deriv3, double_cos_deriv4, double_cos_deriv5};
    return jet5_apply(a, double_cos, derivs);
}

static inline jet5
jet5_tan(jet5 a)
{
    static double (*const derivs[5])(double) = {double_tan_deriv1, double_tan_deriv2, double_tan_deriv3, double_tan_deriv4, double_tan_deriv5};
    return jet5_apply(a, double_tan, derivs);
}

static inline jet5
jet5_asin(jet5 a)
{
    static double (*const derivs[5])(double) = {double_asin_deriv1, double_asin_deriv2, double_asin_deriv3, double_asin_deriv4, double_asin_deriv5};
    return jet5_apply(a, double_asin, derivs);
}

static inline jet5
jet5_acos(jet5 a)
{
    static double (*const derivs[5])(double) = {double_acos_deriv1, double_acos_deriv2, double_acos_deriv3, double_acos_deriv4, double_acos_deriv5};
    return jet5_apply(a, double_acos, derivs);
}

static inline jet5
jet5_atan(jet5 a)
{
    static double (*const derivs[5])(double) = {double_atan_deriv1, double_atan_deriv2, double_atan_deriv3, double_atan_deriv4, double_atan_deriv5};
    return jet5_apply(a, double_atan, derivs);
}

static inline jet5
jet5_sinh(jet5 a)
{
    static double (*const derivs[5])(double) = {double_sinh_deriv1, double_sinh_deriv2, double_sinh_deriv3, double_sinh_deriv4, double_sinh_deriv5};
    return jet5_apply(a, double_sinh, derivs);
}

static inline jet5
jet5_cosh(jet5 a)
{
    static double (*const derivs[5])(double) = {double_cosh_deriv1, double_cosh_deriv2, double_cosh_deriv3, double_cosh_deriv4, double_cosh_deriv5};
    return jet5_apply(a, double_cosh, derivs);
}

static inline jet5
jet5_tanh(jet5 a)
{
    static double (*const derivs[5])(double) = {double_tanh_deriv1, double_tanh_deriv2, double_tanh_deriv3, double_tanh_deriv4, double_tanh_deriv5};
    return jet5_apply(a, double_tanh, derivs);
}

static inline jet5
jet5_asinh(jet5 a)
{
    static double (*const derivs[5])(double) = {double_asinh_deriv1, double_asinh_deriv2, double_asinh_deriv3, double_asinh_deriv4, double_asinh_deriv5};
    return jet5_apply(a, double_asinh, derivs);
}

static inline jet5
jet5_acosh(jet5 a)
{
    static double (*const derivs[5])(double) = {double_acosh_deriv1, double_acosh_deriv2, double_acosh_deriv3, double_acosh_deriv4, double_acosh_deriv5};
    return jet5_apply(a, double_acosh, derivs);
}

static inline jet5
jet5_atanh(jet5 a)
{
    static double (*const derivs[5])(double) = {double_atanh_deriv1, double_atanh_deriv2, double_atanh_deriv3, double_atanh_deriv4, double_atanh_deriv5};
    return jet5_apply(a, double_atanh, derivs);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline jet5
jet5_ipow(jet5 a, int64_t b)
{
    jet5 result = jet5_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return jet5_square(a);
    }
    if (b == -1) {
        return jet5_inv(a);
    }
    if (b < 0) {
        a = jet5_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = jet5_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = jet5_square(a);
        }
    }
    return result;
}

/* Real analytic power: requires a positive real part for log(a).
 * Jets are commutative, so lpow and rpow agree. */
static inline jet5
jet5_rpow(jet5 a, jet5 b)
{
    return jet5_exp(jet5_mul(jet5_log(a), b));
}

static inline jet5
jet5_pow(jet5 a, jet5 b)
{
    return jet5_rpow(a, b);
}

static inline jet5
jet5_lpow(jet5 a, jet5 b)
{
    return jet5_rpow(a, b);
}

static inline jet5
jet5_radd(double a, jet5 b)
{
    return jet5_add(jet5_from_re(a), b);
}

static inline jet5
jet5_ladd(jet5 a, double b)
{
    return jet5_add(a, jet5_from_re(b));
}

static inline jet5
jet5_rsub(double a, jet5 b)
{
    return jet5_sub(jet5_from_re(a), b);
}

static inline jet5
jet5_lsub(jet5 a, double b)
{
    return jet5_sub(a, jet5_from_re(b));
}

static inline jet5
jet5_rmul(double a, jet5 b)
{
    return vec6_scale(b, a);
}

static inline jet5
jet5_lmul(jet5 a, double b)
{
    return vec6_scale(a, b);
}

static inline jet5
jet5_rdiv(double a, jet5 b)
{
    return vec6_scale(jet5_inv(b), a);
}

static inline jet5
jet5_ldiv(jet5 a, double b)
{
    return vec6_div_scalar(a, b);
}

static inline jet5
jet5_rpows(double a, jet5 b)
{
    return jet5_rpow(jet5_from_re(a), b);
}

static inline jet5
jet5_lpows(jet5 a, double b)
{
    return jet5_rpow(a, jet5_from_re(b));
}

static inline jet5
jet5_mix(jet5 a, jet5 b, double t)
{
    jet5 delta = jet5_mul(jet5_inv(a), b);
    return jet5_mul(a, jet5_lpows(delta, t));
}

/* |det lmat(a)|^(1/6): the absolute value of the real part. */
static inline double
jet5_norm(jet5 a)
{
    return fabs(a.data[0]);
}

/* Generic minimal polynomial (t - x)^6, ascending coefficients. */
static inline void
jet5_minpoly(jet5 a, double result[7])
{
    double x = a.data[0];

    result[0] = x*x*x*x*x*x;
    result[1] = -6*x*x*x*x*x;
    result[2] = 15*x*x*x*x;
    result[3] = -20*x*x*x;
    result[4] = 15*x*x;
    result[5] = -6*x;
    result[6] = 1;
}

static inline void
jet5_eivals(jet5 a, complex result[6])
{
    int i;
    for (i = 0; i < 6; i++) {
        result[i] = complex_make((const double[]){a.data[0], 0});
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
jet5_lmat(jet5 a, double result[ALGPY_JET5_N*ALGPY_JET5_N])
{
    jet5 basis = {{0}};
    jet5 product;
    int row, col;
    for (col = 0; col < ALGPY_JET5_N; col++) {
        basis.data[col] = 1;
        product = jet5_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_JET5_N; row++) {
            result[row * ALGPY_JET5_N + col] = product.data[row];
        }
    }
}

static inline void
jet5_rmat(jet5 a, double result[ALGPY_JET5_N*ALGPY_JET5_N])
{
    jet5 basis = {{0}};
    jet5 product;
    int row, col;
    for (row = 0; row < ALGPY_JET5_N; row++) {
        basis.data[row] = 1;
        product = jet5_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_JET5_N; col++) {
            result[row * ALGPY_JET5_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_JET5_H */
