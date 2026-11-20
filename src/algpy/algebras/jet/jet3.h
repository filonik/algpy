/*
 * Jets of order 3: truncated polynomials
 * x + x1*e + ... + x3*e^3 in R[e]/(e^4), stored as a vec4.
 * A function applied to a = x + e gives its Taylor coefficients, the
 * normalized derivatives f^(k)(x)/k!, k = 0..3.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_JET3_H
#define ALGPY_JET3_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_JET3_N 4

typedef vec4 jet3;

static inline jet3
jet3_make(const double coeffs[ALGPY_JET3_N])
{
    jet3 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline jet3
jet3_from_re(double a)
{
    return jet3_make((const double[ALGPY_JET3_N]){a});
}

static inline double
jet3_to_re(jet3 a)
{
    return a.data[0];
}

static inline jet3
jet3_zero(void)
{
    return vec4_zero();
}

static inline jet3
jet3_one(void)
{
    return jet3_from_re(1);
}

static inline int
jet3_is_real(jet3 a)
{
    return vec4_is_real(a);
}

static inline jet3
jet3_re(jet3 a)
{
    return jet3_from_re(a.data[0]);
}

static inline jet3
jet3_im(jet3 a)
{
    jet3 result = a;
    result.data[0] = 0;
    return result;
}

static inline jet3
jet3_pos(jet3 a)
{
    return vec4_pos(a);
}

static inline jet3
jet3_neg(jet3 a)
{
    return vec4_neg(a);
}

static inline jet3
jet3_add(jet3 a, jet3 b)
{
    return vec4_add(a, b);
}

static inline jet3
jet3_sub(jet3 a, jet3 b)
{
    return vec4_sub(a, b);
}

static inline int
jet3_equal(jet3 a, jet3 b)
{
    return vec4_all_eq(a, b);
}

static inline int
jet3_not_equal(jet3 a, jet3 b)
{
    return !jet3_equal(a, b);
}

/* Truncated Cauchy product: coefficient k is the sum of x[i]*y[k-i]. */
static inline jet3
jet3_mul(jet3 a, jet3 b)
{
    const double *x = a.data, *y = b.data;
    return jet3_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] + x[1]*y[1] + x[2]*y[0],
        x[0]*y[3] + x[1]*y[2] + x[2]*y[1] + x[3]*y[0]});
}

static inline jet3
jet3_square(jet3 a)
{
    const double *x = a.data;
    return jet3_make((const double[]){
        x[0]*x[0],
        2*x[0]*x[1],
        2*x[0]*x[2] + x[1]*x[1],
        2*x[0]*x[3] + 2*x[1]*x[2]});
}

/* The automorphism e -> -e: coefficient k changes sign for odd k. */
static inline jet3
jet3_conjugate(jet3 a)
{
    const double *x = a.data;
    return jet3_make((const double[]){x[0], -x[1], x[2], -x[3]});
}

static inline void
jet3_basis(jet3 result[ALGPY_JET3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_JET3_N; i++) {
        for (j = 0; j < ALGPY_JET3_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

/* The single generator e, at coefficient position 1. */
static inline void
jet3_gens(jet3 result[1])
{
    result[0] = jet3_make((const double[ALGPY_JET3_N]){[1] = 1});
}

static inline void
jet3_coeffs(jet3 result[ALGPY_JET3_N * ALGPY_JET3_N])
{
    jet3 basis[ALGPY_JET3_N];
    int i, j;
    jet3_basis(basis);
    for (i = 0; i < ALGPY_JET3_N; i++) {
        for (j = 0; j < ALGPY_JET3_N; j++) {
            result[i * ALGPY_JET3_N + j] = jet3_mul(basis[i], basis[j]);
        }
    }
}

/* Grade projections: grade k is the coefficient of e^k. */
static inline jet3
jet3_grade0(jet3 a)
{
    const double *x = a.data;
    return jet3_make((const double[]){x[0], 0, 0, 0});
}

static inline jet3
jet3_grade1(jet3 a)
{
    const double *x = a.data;
    return jet3_make((const double[]){0, x[1], 0, 0});
}

static inline jet3
jet3_grade2(jet3 a)
{
    const double *x = a.data;
    return jet3_make((const double[]){0, 0, x[2], 0});
}

static inline jet3
jet3_grade3(jet3 a)
{
    const double *x = a.data;
    return jet3_make((const double[]){0, 0, 0, x[3]});
}

/* Reciprocal by the series recurrence y0 = 1/x0 and
 * y_k = -(x1*y_(k-1) + ... + xk*y0)/x0. Requires a nonzero real part. */
static inline jet3
jet3_inv(jet3 a)
{
    const double *x = a.data;
    jet3 result;
    double *y = result.data;
    int k, i;
    y[0] = 1/x[0];
    for (k = 1; k < ALGPY_JET3_N; k++) {
        double sum = 0;
        for (i = 1; i <= k; i++) {
            sum += x[i]*y[k - i];
        }
        y[k] = -sum/x[0];
    }
    return result;
}

static inline jet3
jet3_div(jet3 a, jet3 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return jet3_mul(a, jet3_inv(b));
}

static const double jet3_factorials[4] = {1., 1., 2., 6.};

/*
 * Taylor expansion about the real part x: with a = x + n, n nilpotent,
 * f(a) = sum over k of f^(k)(x)/k! * n^k, evaluated by Horner's rule.
 * derivs[k-1] is the k-th derivative; f must be 3 times
 * differentiable at x. Singular derivatives follow IEEE propagation.
 */
static inline jet3
jet3_apply(jet3 a, double (*f)(double), double (*const derivs[3])(double))
{
    double x = a.data[0];
    jet3 n = a, result = {{0}};
    int k;
    n.data[0] = 0;
    result.data[0] = derivs[2](x)/6.;
    for (k = 2; k >= 1; k--) {
        result = jet3_mul(result, n);
        result.data[0] += derivs[k - 1](x)/jet3_factorials[k];
    }
    result = jet3_mul(result, n);
    result.data[0] += f(x);
    return result;
}

static inline jet3
jet3_exp(jet3 a)
{
    static double (*const derivs[3])(double) = {double_exp_deriv1, double_exp_deriv2, double_exp_deriv3};
    return jet3_apply(a, double_exp, derivs);
}

static inline jet3
jet3_expm1(jet3 a)
{
    static double (*const derivs[3])(double) = {double_expm1_deriv1, double_expm1_deriv2, double_expm1_deriv3};
    return jet3_apply(a, double_expm1, derivs);
}

static inline jet3
jet3_log(jet3 a)
{
    static double (*const derivs[3])(double) = {double_log_deriv1, double_log_deriv2, double_log_deriv3};
    return jet3_apply(a, double_log, derivs);
}

static inline jet3
jet3_log1p(jet3 a)
{
    static double (*const derivs[3])(double) = {double_log1p_deriv1, double_log1p_deriv2, double_log1p_deriv3};
    return jet3_apply(a, double_log1p, derivs);
}

static inline jet3
jet3_sqrt(jet3 a)
{
    static double (*const derivs[3])(double) = {double_sqrt_deriv1, double_sqrt_deriv2, double_sqrt_deriv3};
    return jet3_apply(a, double_sqrt, derivs);
}

static inline jet3
jet3_sin(jet3 a)
{
    static double (*const derivs[3])(double) = {double_sin_deriv1, double_sin_deriv2, double_sin_deriv3};
    return jet3_apply(a, double_sin, derivs);
}

static inline jet3
jet3_cos(jet3 a)
{
    static double (*const derivs[3])(double) = {double_cos_deriv1, double_cos_deriv2, double_cos_deriv3};
    return jet3_apply(a, double_cos, derivs);
}

static inline jet3
jet3_tan(jet3 a)
{
    static double (*const derivs[3])(double) = {double_tan_deriv1, double_tan_deriv2, double_tan_deriv3};
    return jet3_apply(a, double_tan, derivs);
}

static inline jet3
jet3_asin(jet3 a)
{
    static double (*const derivs[3])(double) = {double_asin_deriv1, double_asin_deriv2, double_asin_deriv3};
    return jet3_apply(a, double_asin, derivs);
}

static inline jet3
jet3_acos(jet3 a)
{
    static double (*const derivs[3])(double) = {double_acos_deriv1, double_acos_deriv2, double_acos_deriv3};
    return jet3_apply(a, double_acos, derivs);
}

static inline jet3
jet3_atan(jet3 a)
{
    static double (*const derivs[3])(double) = {double_atan_deriv1, double_atan_deriv2, double_atan_deriv3};
    return jet3_apply(a, double_atan, derivs);
}

static inline jet3
jet3_sinh(jet3 a)
{
    static double (*const derivs[3])(double) = {double_sinh_deriv1, double_sinh_deriv2, double_sinh_deriv3};
    return jet3_apply(a, double_sinh, derivs);
}

static inline jet3
jet3_cosh(jet3 a)
{
    static double (*const derivs[3])(double) = {double_cosh_deriv1, double_cosh_deriv2, double_cosh_deriv3};
    return jet3_apply(a, double_cosh, derivs);
}

static inline jet3
jet3_tanh(jet3 a)
{
    static double (*const derivs[3])(double) = {double_tanh_deriv1, double_tanh_deriv2, double_tanh_deriv3};
    return jet3_apply(a, double_tanh, derivs);
}

static inline jet3
jet3_asinh(jet3 a)
{
    static double (*const derivs[3])(double) = {double_asinh_deriv1, double_asinh_deriv2, double_asinh_deriv3};
    return jet3_apply(a, double_asinh, derivs);
}

static inline jet3
jet3_acosh(jet3 a)
{
    static double (*const derivs[3])(double) = {double_acosh_deriv1, double_acosh_deriv2, double_acosh_deriv3};
    return jet3_apply(a, double_acosh, derivs);
}

static inline jet3
jet3_atanh(jet3 a)
{
    static double (*const derivs[3])(double) = {double_atanh_deriv1, double_atanh_deriv2, double_atanh_deriv3};
    return jet3_apply(a, double_atanh, derivs);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline jet3
jet3_ipow(jet3 a, int64_t b)
{
    jet3 result = jet3_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return jet3_square(a);
    }
    if (b == -1) {
        return jet3_inv(a);
    }
    if (b < 0) {
        a = jet3_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = jet3_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = jet3_square(a);
        }
    }
    return result;
}

/* Real analytic power: requires a positive real part for log(a).
 * Jets are commutative, so lpow and rpow agree. */
static inline jet3
jet3_rpow(jet3 a, jet3 b)
{
    return jet3_exp(jet3_mul(jet3_log(a), b));
}

static inline jet3
jet3_pow(jet3 a, jet3 b)
{
    return jet3_rpow(a, b);
}

static inline jet3
jet3_lpow(jet3 a, jet3 b)
{
    return jet3_rpow(a, b);
}

static inline jet3
jet3_radd(double a, jet3 b)
{
    return jet3_add(jet3_from_re(a), b);
}

static inline jet3
jet3_ladd(jet3 a, double b)
{
    return jet3_add(a, jet3_from_re(b));
}

static inline jet3
jet3_rsub(double a, jet3 b)
{
    return jet3_sub(jet3_from_re(a), b);
}

static inline jet3
jet3_lsub(jet3 a, double b)
{
    return jet3_sub(a, jet3_from_re(b));
}

static inline jet3
jet3_rmul(double a, jet3 b)
{
    return vec4_scale(b, a);
}

static inline jet3
jet3_lmul(jet3 a, double b)
{
    return vec4_scale(a, b);
}

static inline jet3
jet3_rdiv(double a, jet3 b)
{
    return vec4_scale(jet3_inv(b), a);
}

static inline jet3
jet3_ldiv(jet3 a, double b)
{
    return vec4_div_scalar(a, b);
}

static inline jet3
jet3_rpows(double a, jet3 b)
{
    return jet3_rpow(jet3_from_re(a), b);
}

static inline jet3
jet3_lpows(jet3 a, double b)
{
    return jet3_rpow(a, jet3_from_re(b));
}

static inline jet3
jet3_mix(jet3 a, jet3 b, double t)
{
    jet3 delta = jet3_mul(jet3_inv(a), b);
    return jet3_mul(a, jet3_lpows(delta, t));
}

/* |det lmat(a)|^(1/4): the absolute value of the real part. */
static inline double
jet3_norm(jet3 a)
{
    return fabs(a.data[0]);
}

/* Generic minimal polynomial (t - x)^4, ascending coefficients. */
static inline void
jet3_minpoly(jet3 a, double result[5])
{
    double x = a.data[0];

    result[0] = x*x*x*x;
    result[1] = -4*x*x*x;
    result[2] = 6*x*x;
    result[3] = -4*x;
    result[4] = 1;
}

static inline void
jet3_eivals(jet3 a, complex result[4])
{
    int i;
    for (i = 0; i < 4; i++) {
        result[i] = complex_make((const double[]){a.data[0], 0});
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
jet3_lmat(jet3 a, double result[ALGPY_JET3_N*ALGPY_JET3_N])
{
    jet3 basis = {{0}};
    jet3 product;
    int row, col;
    for (col = 0; col < ALGPY_JET3_N; col++) {
        basis.data[col] = 1;
        product = jet3_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_JET3_N; row++) {
            result[row * ALGPY_JET3_N + col] = product.data[row];
        }
    }
}

static inline void
jet3_rmat(jet3 a, double result[ALGPY_JET3_N*ALGPY_JET3_N])
{
    jet3 basis = {{0}};
    jet3 product;
    int row, col;
    for (row = 0; row < ALGPY_JET3_N; row++) {
        basis.data[row] = 1;
        product = jet3_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_JET3_N; col++) {
            result[row * ALGPY_JET3_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_JET3_H */
