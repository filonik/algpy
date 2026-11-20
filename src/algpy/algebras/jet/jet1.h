/*
 * Jets of order 1: truncated polynomials
 * x + x1*e + ... + x1*e^1 in R[e]/(e^2), stored as a vec2.
 * A function applied to a = x + e gives its Taylor coefficients, the
 * normalized derivatives f^(k)(x)/k!, k = 0..1.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_JET1_H
#define ALGPY_JET1_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../vectorspaces/vec2.h"

#define ALGPY_JET1_N 2

typedef vec2 jet1;

static inline jet1
jet1_make(const double coeffs[ALGPY_JET1_N])
{
    jet1 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline jet1
jet1_from_re(double a)
{
    return jet1_make((const double[ALGPY_JET1_N]){a});
}

static inline double
jet1_to_re(jet1 a)
{
    return a.data[0];
}

static inline jet1
jet1_zero(void)
{
    return vec2_zero();
}

static inline jet1
jet1_one(void)
{
    return jet1_from_re(1);
}

static inline int
jet1_is_real(jet1 a)
{
    return vec2_is_real(a);
}

static inline jet1
jet1_re(jet1 a)
{
    return jet1_from_re(a.data[0]);
}

static inline jet1
jet1_im(jet1 a)
{
    jet1 result = a;
    result.data[0] = 0;
    return result;
}

static inline jet1
jet1_pos(jet1 a)
{
    return vec2_pos(a);
}

static inline jet1
jet1_neg(jet1 a)
{
    return vec2_neg(a);
}

static inline jet1
jet1_add(jet1 a, jet1 b)
{
    return vec2_add(a, b);
}

static inline jet1
jet1_sub(jet1 a, jet1 b)
{
    return vec2_sub(a, b);
}

static inline int
jet1_equal(jet1 a, jet1 b)
{
    return vec2_all_eq(a, b);
}

static inline int
jet1_not_equal(jet1 a, jet1 b)
{
    return !jet1_equal(a, b);
}

/* Truncated Cauchy product: coefficient k is the sum of x[i]*y[k-i]. */
static inline jet1
jet1_mul(jet1 a, jet1 b)
{
    const double *x = a.data, *y = b.data;
    return jet1_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0]});
}

static inline jet1
jet1_square(jet1 a)
{
    const double *x = a.data;
    return jet1_make((const double[]){
        x[0]*x[0],
        2*x[0]*x[1]});
}

/* The automorphism e -> -e: coefficient k changes sign for odd k. */
static inline jet1
jet1_conjugate(jet1 a)
{
    const double *x = a.data;
    return jet1_make((const double[]){x[0], -x[1]});
}

static inline void
jet1_basis(jet1 result[ALGPY_JET1_N])
{
    int i, j;
    for (i = 0; i < ALGPY_JET1_N; i++) {
        for (j = 0; j < ALGPY_JET1_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

/* The single generator e, at coefficient position 1. */
static inline void
jet1_gens(jet1 result[1])
{
    result[0] = jet1_make((const double[ALGPY_JET1_N]){[1] = 1});
}

static inline void
jet1_coeffs(jet1 result[ALGPY_JET1_N * ALGPY_JET1_N])
{
    jet1 basis[ALGPY_JET1_N];
    int i, j;
    jet1_basis(basis);
    for (i = 0; i < ALGPY_JET1_N; i++) {
        for (j = 0; j < ALGPY_JET1_N; j++) {
            result[i * ALGPY_JET1_N + j] = jet1_mul(basis[i], basis[j]);
        }
    }
}

/* Grade projections: grade k is the coefficient of e^k. */
static inline jet1
jet1_grade0(jet1 a)
{
    const double *x = a.data;
    return jet1_make((const double[]){x[0], 0});
}

static inline jet1
jet1_grade1(jet1 a)
{
    const double *x = a.data;
    return jet1_make((const double[]){0, x[1]});
}

/* Reciprocal by the series recurrence y0 = 1/x0 and
 * y_k = -(x1*y_(k-1) + ... + xk*y0)/x0. Requires a nonzero real part. */
static inline jet1
jet1_inv(jet1 a)
{
    const double *x = a.data;
    jet1 result;
    double *y = result.data;
    int k, i;
    y[0] = 1/x[0];
    for (k = 1; k < ALGPY_JET1_N; k++) {
        double sum = 0;
        for (i = 1; i <= k; i++) {
            sum += x[i]*y[k - i];
        }
        y[k] = -sum/x[0];
    }
    return result;
}

static inline jet1
jet1_div(jet1 a, jet1 b)
{
    if (vec2_is_real(b)) {
        return vec2_div_scalar(a, b.data[0]);
    }
    return jet1_mul(a, jet1_inv(b));
}

static const double jet1_factorials[2] = {1., 1.};

/*
 * Taylor expansion about the real part x: with a = x + n, n nilpotent,
 * f(a) = sum over k of f^(k)(x)/k! * n^k, evaluated by Horner's rule.
 * derivs[k-1] is the k-th derivative; f must be 1 times
 * differentiable at x. Singular derivatives follow IEEE propagation.
 */
static inline jet1
jet1_apply(jet1 a, double (*f)(double), double (*const derivs[1])(double))
{
    double x = a.data[0];
    jet1 n = a, result = {{0}};
    int k;
    n.data[0] = 0;
    result.data[0] = derivs[0](x)/1.;
    for (k = 0; k >= 1; k--) {
        result = jet1_mul(result, n);
        result.data[0] += derivs[k - 1](x)/jet1_factorials[k];
    }
    result = jet1_mul(result, n);
    result.data[0] += f(x);
    return result;
}

static inline jet1
jet1_exp(jet1 a)
{
    static double (*const derivs[1])(double) = {double_exp_deriv1};
    return jet1_apply(a, double_exp, derivs);
}

static inline jet1
jet1_expm1(jet1 a)
{
    static double (*const derivs[1])(double) = {double_expm1_deriv1};
    return jet1_apply(a, double_expm1, derivs);
}

static inline jet1
jet1_log(jet1 a)
{
    static double (*const derivs[1])(double) = {double_log_deriv1};
    return jet1_apply(a, double_log, derivs);
}

static inline jet1
jet1_log1p(jet1 a)
{
    static double (*const derivs[1])(double) = {double_log1p_deriv1};
    return jet1_apply(a, double_log1p, derivs);
}

static inline jet1
jet1_sqrt(jet1 a)
{
    static double (*const derivs[1])(double) = {double_sqrt_deriv1};
    return jet1_apply(a, double_sqrt, derivs);
}

static inline jet1
jet1_sin(jet1 a)
{
    static double (*const derivs[1])(double) = {double_sin_deriv1};
    return jet1_apply(a, double_sin, derivs);
}

static inline jet1
jet1_cos(jet1 a)
{
    static double (*const derivs[1])(double) = {double_cos_deriv1};
    return jet1_apply(a, double_cos, derivs);
}

static inline jet1
jet1_tan(jet1 a)
{
    static double (*const derivs[1])(double) = {double_tan_deriv1};
    return jet1_apply(a, double_tan, derivs);
}

static inline jet1
jet1_asin(jet1 a)
{
    static double (*const derivs[1])(double) = {double_asin_deriv1};
    return jet1_apply(a, double_asin, derivs);
}

static inline jet1
jet1_acos(jet1 a)
{
    static double (*const derivs[1])(double) = {double_acos_deriv1};
    return jet1_apply(a, double_acos, derivs);
}

static inline jet1
jet1_atan(jet1 a)
{
    static double (*const derivs[1])(double) = {double_atan_deriv1};
    return jet1_apply(a, double_atan, derivs);
}

static inline jet1
jet1_sinh(jet1 a)
{
    static double (*const derivs[1])(double) = {double_sinh_deriv1};
    return jet1_apply(a, double_sinh, derivs);
}

static inline jet1
jet1_cosh(jet1 a)
{
    static double (*const derivs[1])(double) = {double_cosh_deriv1};
    return jet1_apply(a, double_cosh, derivs);
}

static inline jet1
jet1_tanh(jet1 a)
{
    static double (*const derivs[1])(double) = {double_tanh_deriv1};
    return jet1_apply(a, double_tanh, derivs);
}

static inline jet1
jet1_asinh(jet1 a)
{
    static double (*const derivs[1])(double) = {double_asinh_deriv1};
    return jet1_apply(a, double_asinh, derivs);
}

static inline jet1
jet1_acosh(jet1 a)
{
    static double (*const derivs[1])(double) = {double_acosh_deriv1};
    return jet1_apply(a, double_acosh, derivs);
}

static inline jet1
jet1_atanh(jet1 a)
{
    static double (*const derivs[1])(double) = {double_atanh_deriv1};
    return jet1_apply(a, double_atanh, derivs);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline jet1
jet1_ipow(jet1 a, int64_t b)
{
    jet1 result = jet1_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return jet1_square(a);
    }
    if (b == -1) {
        return jet1_inv(a);
    }
    if (b < 0) {
        a = jet1_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = jet1_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = jet1_square(a);
        }
    }
    return result;
}

/* Real analytic power: requires a positive real part for log(a).
 * Jets are commutative, so lpow and rpow agree. */
static inline jet1
jet1_rpow(jet1 a, jet1 b)
{
    return jet1_exp(jet1_mul(jet1_log(a), b));
}

static inline jet1
jet1_pow(jet1 a, jet1 b)
{
    return jet1_rpow(a, b);
}

static inline jet1
jet1_lpow(jet1 a, jet1 b)
{
    return jet1_rpow(a, b);
}

static inline jet1
jet1_radd(double a, jet1 b)
{
    return jet1_add(jet1_from_re(a), b);
}

static inline jet1
jet1_ladd(jet1 a, double b)
{
    return jet1_add(a, jet1_from_re(b));
}

static inline jet1
jet1_rsub(double a, jet1 b)
{
    return jet1_sub(jet1_from_re(a), b);
}

static inline jet1
jet1_lsub(jet1 a, double b)
{
    return jet1_sub(a, jet1_from_re(b));
}

static inline jet1
jet1_rmul(double a, jet1 b)
{
    return vec2_scale(b, a);
}

static inline jet1
jet1_lmul(jet1 a, double b)
{
    return vec2_scale(a, b);
}

static inline jet1
jet1_rdiv(double a, jet1 b)
{
    return vec2_scale(jet1_inv(b), a);
}

static inline jet1
jet1_ldiv(jet1 a, double b)
{
    return vec2_div_scalar(a, b);
}

static inline jet1
jet1_rpows(double a, jet1 b)
{
    return jet1_rpow(jet1_from_re(a), b);
}

static inline jet1
jet1_lpows(jet1 a, double b)
{
    return jet1_rpow(a, jet1_from_re(b));
}

static inline jet1
jet1_mix(jet1 a, jet1 b, double t)
{
    jet1 delta = jet1_mul(jet1_inv(a), b);
    return jet1_mul(a, jet1_lpows(delta, t));
}

/* |det lmat(a)|^(1/2): the absolute value of the real part. */
static inline double
jet1_norm(jet1 a)
{
    return fabs(a.data[0]);
}

/* Generic minimal polynomial (t - x)^2, ascending coefficients. */
static inline void
jet1_minpoly(jet1 a, double result[3])
{
    double x = a.data[0];

    result[0] = x*x;
    result[1] = -2*x;
    result[2] = 1;
}

static inline void
jet1_eivals(jet1 a, complex result[2])
{
    int i;
    for (i = 0; i < 2; i++) {
        result[i] = complex_make((const double[]){a.data[0], 0});
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
jet1_lmat(jet1 a, double result[ALGPY_JET1_N*ALGPY_JET1_N])
{
    jet1 basis = {{0}};
    jet1 product;
    int row, col;
    for (col = 0; col < ALGPY_JET1_N; col++) {
        basis.data[col] = 1;
        product = jet1_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_JET1_N; row++) {
            result[row * ALGPY_JET1_N + col] = product.data[row];
        }
    }
}

static inline void
jet1_rmat(jet1 a, double result[ALGPY_JET1_N*ALGPY_JET1_N])
{
    jet1 basis = {{0}};
    jet1 product;
    int row, col;
    for (row = 0; row < ALGPY_JET1_N; row++) {
        basis.data[row] = 1;
        product = jet1_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_JET1_N; col++) {
            result[row * ALGPY_JET1_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_JET1_H */
