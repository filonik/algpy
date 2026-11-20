/*
 * Jets of order 4: truncated polynomials
 * x + x1*e + ... + x4*e^4 in R[e]/(e^5), stored as a vec5.
 * A function applied to a = x + e gives its Taylor coefficients, the
 * normalized derivatives f^(k)(x)/k!, k = 0..4.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_JET4_H
#define ALGPY_JET4_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../vectorspaces/vec5.h"

#define ALGPY_JET4_N 5

typedef vec5 jet4;

static inline jet4
jet4_make(const double coeffs[ALGPY_JET4_N])
{
    jet4 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline jet4
jet4_from_re(double a)
{
    return jet4_make((const double[ALGPY_JET4_N]){a});
}

static inline double
jet4_to_re(jet4 a)
{
    return a.data[0];
}

static inline jet4
jet4_zero(void)
{
    return vec5_zero();
}

static inline jet4
jet4_one(void)
{
    return jet4_from_re(1);
}

static inline int
jet4_is_real(jet4 a)
{
    return vec5_is_real(a);
}

static inline jet4
jet4_re(jet4 a)
{
    return jet4_from_re(a.data[0]);
}

static inline jet4
jet4_im(jet4 a)
{
    jet4 result = a;
    result.data[0] = 0;
    return result;
}

static inline jet4
jet4_pos(jet4 a)
{
    return vec5_pos(a);
}

static inline jet4
jet4_neg(jet4 a)
{
    return vec5_neg(a);
}

static inline jet4
jet4_add(jet4 a, jet4 b)
{
    return vec5_add(a, b);
}

static inline jet4
jet4_sub(jet4 a, jet4 b)
{
    return vec5_sub(a, b);
}

static inline int
jet4_equal(jet4 a, jet4 b)
{
    return vec5_all_eq(a, b);
}

static inline int
jet4_not_equal(jet4 a, jet4 b)
{
    return !jet4_equal(a, b);
}

/* Truncated Cauchy product: coefficient k is the sum of x[i]*y[k-i]. */
static inline jet4
jet4_mul(jet4 a, jet4 b)
{
    const double *x = a.data, *y = b.data;
    return jet4_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] + x[1]*y[1] + x[2]*y[0],
        x[0]*y[3] + x[1]*y[2] + x[2]*y[1] + x[3]*y[0],
        x[0]*y[4] + x[1]*y[3] + x[2]*y[2] + x[3]*y[1] + x[4]*y[0]});
}

static inline jet4
jet4_square(jet4 a)
{
    const double *x = a.data;
    return jet4_make((const double[]){
        x[0]*x[0],
        2*x[0]*x[1],
        2*x[0]*x[2] + x[1]*x[1],
        2*x[0]*x[3] + 2*x[1]*x[2],
        2*x[0]*x[4] + 2*x[1]*x[3] + x[2]*x[2]});
}

/* The automorphism e -> -e: coefficient k changes sign for odd k. */
static inline jet4
jet4_conjugate(jet4 a)
{
    const double *x = a.data;
    return jet4_make((const double[]){x[0], -x[1], x[2], -x[3], x[4]});
}

static inline void
jet4_basis(jet4 result[ALGPY_JET4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_JET4_N; i++) {
        for (j = 0; j < ALGPY_JET4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

/* The single generator e, at coefficient position 1. */
static inline void
jet4_gens(jet4 result[1])
{
    result[0] = jet4_make((const double[ALGPY_JET4_N]){[1] = 1});
}

static inline void
jet4_coeffs(jet4 result[ALGPY_JET4_N * ALGPY_JET4_N])
{
    jet4 basis[ALGPY_JET4_N];
    int i, j;
    jet4_basis(basis);
    for (i = 0; i < ALGPY_JET4_N; i++) {
        for (j = 0; j < ALGPY_JET4_N; j++) {
            result[i * ALGPY_JET4_N + j] = jet4_mul(basis[i], basis[j]);
        }
    }
}

/* Grade projections: grade k is the coefficient of e^k. */
static inline jet4
jet4_grade0(jet4 a)
{
    const double *x = a.data;
    return jet4_make((const double[]){x[0], 0, 0, 0, 0});
}

static inline jet4
jet4_grade1(jet4 a)
{
    const double *x = a.data;
    return jet4_make((const double[]){0, x[1], 0, 0, 0});
}

static inline jet4
jet4_grade2(jet4 a)
{
    const double *x = a.data;
    return jet4_make((const double[]){0, 0, x[2], 0, 0});
}

static inline jet4
jet4_grade3(jet4 a)
{
    const double *x = a.data;
    return jet4_make((const double[]){0, 0, 0, x[3], 0});
}

static inline jet4
jet4_grade4(jet4 a)
{
    const double *x = a.data;
    return jet4_make((const double[]){0, 0, 0, 0, x[4]});
}

/* Reciprocal by the series recurrence y0 = 1/x0 and
 * y_k = -(x1*y_(k-1) + ... + xk*y0)/x0. Requires a nonzero real part. */
static inline jet4
jet4_inv(jet4 a)
{
    const double *x = a.data;
    jet4 result;
    double *y = result.data;
    int k, i;
    y[0] = 1/x[0];
    for (k = 1; k < ALGPY_JET4_N; k++) {
        double sum = 0;
        for (i = 1; i <= k; i++) {
            sum += x[i]*y[k - i];
        }
        y[k] = -sum/x[0];
    }
    return result;
}

static inline jet4
jet4_div(jet4 a, jet4 b)
{
    if (vec5_is_real(b)) {
        return vec5_div_scalar(a, b.data[0]);
    }
    return jet4_mul(a, jet4_inv(b));
}

static const double jet4_factorials[5] = {1., 1., 2., 6., 24.};

/*
 * Taylor expansion about the real part x: with a = x + n, n nilpotent,
 * f(a) = sum over k of f^(k)(x)/k! * n^k, evaluated by Horner's rule.
 * derivs[k-1] is the k-th derivative; f must be 4 times
 * differentiable at x. Singular derivatives follow IEEE propagation.
 */
static inline jet4
jet4_apply(jet4 a, double (*f)(double), double (*const derivs[4])(double))
{
    double x = a.data[0];
    jet4 n = a, result = {{0}};
    int k;
    n.data[0] = 0;
    result.data[0] = derivs[3](x)/24.;
    for (k = 3; k >= 1; k--) {
        result = jet4_mul(result, n);
        result.data[0] += derivs[k - 1](x)/jet4_factorials[k];
    }
    result = jet4_mul(result, n);
    result.data[0] += f(x);
    return result;
}

static inline jet4
jet4_exp(jet4 a)
{
    static double (*const derivs[4])(double) = {double_exp_deriv1, double_exp_deriv2, double_exp_deriv3, double_exp_deriv4};
    return jet4_apply(a, double_exp, derivs);
}

static inline jet4
jet4_expm1(jet4 a)
{
    static double (*const derivs[4])(double) = {double_expm1_deriv1, double_expm1_deriv2, double_expm1_deriv3, double_expm1_deriv4};
    return jet4_apply(a, double_expm1, derivs);
}

static inline jet4
jet4_log(jet4 a)
{
    static double (*const derivs[4])(double) = {double_log_deriv1, double_log_deriv2, double_log_deriv3, double_log_deriv4};
    return jet4_apply(a, double_log, derivs);
}

static inline jet4
jet4_log1p(jet4 a)
{
    static double (*const derivs[4])(double) = {double_log1p_deriv1, double_log1p_deriv2, double_log1p_deriv3, double_log1p_deriv4};
    return jet4_apply(a, double_log1p, derivs);
}

static inline jet4
jet4_sqrt(jet4 a)
{
    static double (*const derivs[4])(double) = {double_sqrt_deriv1, double_sqrt_deriv2, double_sqrt_deriv3, double_sqrt_deriv4};
    return jet4_apply(a, double_sqrt, derivs);
}

static inline jet4
jet4_sin(jet4 a)
{
    static double (*const derivs[4])(double) = {double_sin_deriv1, double_sin_deriv2, double_sin_deriv3, double_sin_deriv4};
    return jet4_apply(a, double_sin, derivs);
}

static inline jet4
jet4_cos(jet4 a)
{
    static double (*const derivs[4])(double) = {double_cos_deriv1, double_cos_deriv2, double_cos_deriv3, double_cos_deriv4};
    return jet4_apply(a, double_cos, derivs);
}

static inline jet4
jet4_tan(jet4 a)
{
    static double (*const derivs[4])(double) = {double_tan_deriv1, double_tan_deriv2, double_tan_deriv3, double_tan_deriv4};
    return jet4_apply(a, double_tan, derivs);
}

static inline jet4
jet4_asin(jet4 a)
{
    static double (*const derivs[4])(double) = {double_asin_deriv1, double_asin_deriv2, double_asin_deriv3, double_asin_deriv4};
    return jet4_apply(a, double_asin, derivs);
}

static inline jet4
jet4_acos(jet4 a)
{
    static double (*const derivs[4])(double) = {double_acos_deriv1, double_acos_deriv2, double_acos_deriv3, double_acos_deriv4};
    return jet4_apply(a, double_acos, derivs);
}

static inline jet4
jet4_atan(jet4 a)
{
    static double (*const derivs[4])(double) = {double_atan_deriv1, double_atan_deriv2, double_atan_deriv3, double_atan_deriv4};
    return jet4_apply(a, double_atan, derivs);
}

static inline jet4
jet4_sinh(jet4 a)
{
    static double (*const derivs[4])(double) = {double_sinh_deriv1, double_sinh_deriv2, double_sinh_deriv3, double_sinh_deriv4};
    return jet4_apply(a, double_sinh, derivs);
}

static inline jet4
jet4_cosh(jet4 a)
{
    static double (*const derivs[4])(double) = {double_cosh_deriv1, double_cosh_deriv2, double_cosh_deriv3, double_cosh_deriv4};
    return jet4_apply(a, double_cosh, derivs);
}

static inline jet4
jet4_tanh(jet4 a)
{
    static double (*const derivs[4])(double) = {double_tanh_deriv1, double_tanh_deriv2, double_tanh_deriv3, double_tanh_deriv4};
    return jet4_apply(a, double_tanh, derivs);
}

static inline jet4
jet4_asinh(jet4 a)
{
    static double (*const derivs[4])(double) = {double_asinh_deriv1, double_asinh_deriv2, double_asinh_deriv3, double_asinh_deriv4};
    return jet4_apply(a, double_asinh, derivs);
}

static inline jet4
jet4_acosh(jet4 a)
{
    static double (*const derivs[4])(double) = {double_acosh_deriv1, double_acosh_deriv2, double_acosh_deriv3, double_acosh_deriv4};
    return jet4_apply(a, double_acosh, derivs);
}

static inline jet4
jet4_atanh(jet4 a)
{
    static double (*const derivs[4])(double) = {double_atanh_deriv1, double_atanh_deriv2, double_atanh_deriv3, double_atanh_deriv4};
    return jet4_apply(a, double_atanh, derivs);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline jet4
jet4_ipow(jet4 a, int64_t b)
{
    jet4 result = jet4_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return jet4_square(a);
    }
    if (b == -1) {
        return jet4_inv(a);
    }
    if (b < 0) {
        a = jet4_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = jet4_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = jet4_square(a);
        }
    }
    return result;
}

/* Real analytic power: requires a positive real part for log(a).
 * Jets are commutative, so lpow and rpow agree. */
static inline jet4
jet4_rpow(jet4 a, jet4 b)
{
    return jet4_exp(jet4_mul(jet4_log(a), b));
}

static inline jet4
jet4_pow(jet4 a, jet4 b)
{
    return jet4_rpow(a, b);
}

static inline jet4
jet4_lpow(jet4 a, jet4 b)
{
    return jet4_rpow(a, b);
}

static inline jet4
jet4_radd(double a, jet4 b)
{
    return jet4_add(jet4_from_re(a), b);
}

static inline jet4
jet4_ladd(jet4 a, double b)
{
    return jet4_add(a, jet4_from_re(b));
}

static inline jet4
jet4_rsub(double a, jet4 b)
{
    return jet4_sub(jet4_from_re(a), b);
}

static inline jet4
jet4_lsub(jet4 a, double b)
{
    return jet4_sub(a, jet4_from_re(b));
}

static inline jet4
jet4_rmul(double a, jet4 b)
{
    return vec5_scale(b, a);
}

static inline jet4
jet4_lmul(jet4 a, double b)
{
    return vec5_scale(a, b);
}

static inline jet4
jet4_rdiv(double a, jet4 b)
{
    return vec5_scale(jet4_inv(b), a);
}

static inline jet4
jet4_ldiv(jet4 a, double b)
{
    return vec5_div_scalar(a, b);
}

static inline jet4
jet4_rpows(double a, jet4 b)
{
    return jet4_rpow(jet4_from_re(a), b);
}

static inline jet4
jet4_lpows(jet4 a, double b)
{
    return jet4_rpow(a, jet4_from_re(b));
}

static inline jet4
jet4_mix(jet4 a, jet4 b, double t)
{
    jet4 delta = jet4_mul(jet4_inv(a), b);
    return jet4_mul(a, jet4_lpows(delta, t));
}

/* |det lmat(a)|^(1/5): the absolute value of the real part. */
static inline double
jet4_norm(jet4 a)
{
    return fabs(a.data[0]);
}

/* Generic minimal polynomial (t - x)^5, ascending coefficients. */
static inline void
jet4_minpoly(jet4 a, double result[6])
{
    double x = a.data[0];

    result[0] = -x*x*x*x*x;
    result[1] = 5*x*x*x*x;
    result[2] = -10*x*x*x;
    result[3] = 10*x*x;
    result[4] = -5*x;
    result[5] = 1;
}

static inline void
jet4_eivals(jet4 a, complex result[5])
{
    int i;
    for (i = 0; i < 5; i++) {
        result[i] = complex_make((const double[]){a.data[0], 0});
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
jet4_lmat(jet4 a, double result[ALGPY_JET4_N*ALGPY_JET4_N])
{
    jet4 basis = {{0}};
    jet4 product;
    int row, col;
    for (col = 0; col < ALGPY_JET4_N; col++) {
        basis.data[col] = 1;
        product = jet4_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_JET4_N; row++) {
            result[row * ALGPY_JET4_N + col] = product.data[row];
        }
    }
}

static inline void
jet4_rmat(jet4 a, double result[ALGPY_JET4_N*ALGPY_JET4_N])
{
    jet4 basis = {{0}};
    jet4 product;
    int row, col;
    for (row = 0; row < ALGPY_JET4_N; row++) {
        basis.data[row] = 1;
        product = jet4_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_JET4_N; col++) {
            result[row * ALGPY_JET4_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_JET4_H */
