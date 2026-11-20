/*
 * Jets of order 2: truncated polynomials
 * x + x1*e + ... + x2*e^2 in R[e]/(e^3), stored as a vec3.
 * A function applied to a = x + e gives its Taylor coefficients, the
 * normalized derivatives f^(k)(x)/k!, k = 0..2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_JET2_H
#define ALGPY_JET2_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../vectorspaces/vec3.h"

#define ALGPY_JET2_N 3

typedef vec3 jet2;

static inline jet2
jet2_make(const double coeffs[ALGPY_JET2_N])
{
    jet2 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline jet2
jet2_from_re(double a)
{
    return jet2_make((const double[ALGPY_JET2_N]){a});
}

static inline double
jet2_to_re(jet2 a)
{
    return a.data[0];
}

static inline jet2
jet2_zero(void)
{
    return vec3_zero();
}

static inline jet2
jet2_one(void)
{
    return jet2_from_re(1);
}

static inline int
jet2_is_real(jet2 a)
{
    return vec3_is_real(a);
}

static inline jet2
jet2_re(jet2 a)
{
    return jet2_from_re(a.data[0]);
}

static inline jet2
jet2_im(jet2 a)
{
    jet2 result = a;
    result.data[0] = 0;
    return result;
}

static inline jet2
jet2_pos(jet2 a)
{
    return vec3_pos(a);
}

static inline jet2
jet2_neg(jet2 a)
{
    return vec3_neg(a);
}

static inline jet2
jet2_add(jet2 a, jet2 b)
{
    return vec3_add(a, b);
}

static inline jet2
jet2_sub(jet2 a, jet2 b)
{
    return vec3_sub(a, b);
}

static inline int
jet2_equal(jet2 a, jet2 b)
{
    return vec3_all_eq(a, b);
}

static inline int
jet2_not_equal(jet2 a, jet2 b)
{
    return !jet2_equal(a, b);
}

/* Truncated Cauchy product: coefficient k is the sum of x[i]*y[k-i]. */
static inline jet2
jet2_mul(jet2 a, jet2 b)
{
    const double *x = a.data, *y = b.data;
    return jet2_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] + x[1]*y[1] + x[2]*y[0]});
}

static inline jet2
jet2_square(jet2 a)
{
    const double *x = a.data;
    return jet2_make((const double[]){
        x[0]*x[0],
        2*x[0]*x[1],
        2*x[0]*x[2] + x[1]*x[1]});
}

/* The automorphism e -> -e: coefficient k changes sign for odd k. */
static inline jet2
jet2_conjugate(jet2 a)
{
    const double *x = a.data;
    return jet2_make((const double[]){x[0], -x[1], x[2]});
}

static inline void
jet2_basis(jet2 result[ALGPY_JET2_N])
{
    int i, j;
    for (i = 0; i < ALGPY_JET2_N; i++) {
        for (j = 0; j < ALGPY_JET2_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

/* The single generator e, at coefficient position 1. */
static inline void
jet2_gens(jet2 result[1])
{
    result[0] = jet2_make((const double[ALGPY_JET2_N]){[1] = 1});
}

static inline void
jet2_coeffs(jet2 result[ALGPY_JET2_N * ALGPY_JET2_N])
{
    jet2 basis[ALGPY_JET2_N];
    int i, j;
    jet2_basis(basis);
    for (i = 0; i < ALGPY_JET2_N; i++) {
        for (j = 0; j < ALGPY_JET2_N; j++) {
            result[i * ALGPY_JET2_N + j] = jet2_mul(basis[i], basis[j]);
        }
    }
}

/* Grade projections: grade k is the coefficient of e^k. */
static inline jet2
jet2_grade0(jet2 a)
{
    const double *x = a.data;
    return jet2_make((const double[]){x[0], 0, 0});
}

static inline jet2
jet2_grade1(jet2 a)
{
    const double *x = a.data;
    return jet2_make((const double[]){0, x[1], 0});
}

static inline jet2
jet2_grade2(jet2 a)
{
    const double *x = a.data;
    return jet2_make((const double[]){0, 0, x[2]});
}

/* Reciprocal by the series recurrence y0 = 1/x0 and
 * y_k = -(x1*y_(k-1) + ... + xk*y0)/x0. Requires a nonzero real part. */
static inline jet2
jet2_inv(jet2 a)
{
    const double *x = a.data;
    jet2 result;
    double *y = result.data;
    int k, i;
    y[0] = 1/x[0];
    for (k = 1; k < ALGPY_JET2_N; k++) {
        double sum = 0;
        for (i = 1; i <= k; i++) {
            sum += x[i]*y[k - i];
        }
        y[k] = -sum/x[0];
    }
    return result;
}

static inline jet2
jet2_div(jet2 a, jet2 b)
{
    if (vec3_is_real(b)) {
        return vec3_div_scalar(a, b.data[0]);
    }
    return jet2_mul(a, jet2_inv(b));
}

static const double jet2_factorials[3] = {1., 1., 2.};

/*
 * Taylor expansion about the real part x: with a = x + n, n nilpotent,
 * f(a) = sum over k of f^(k)(x)/k! * n^k, evaluated by Horner's rule.
 * derivs[k-1] is the k-th derivative; f must be 2 times
 * differentiable at x. Singular derivatives follow IEEE propagation.
 */
static inline jet2
jet2_apply(jet2 a, double (*f)(double), double (*const derivs[2])(double))
{
    double x = a.data[0];
    jet2 n = a, result = {{0}};
    int k;
    n.data[0] = 0;
    result.data[0] = derivs[1](x)/2.;
    for (k = 1; k >= 1; k--) {
        result = jet2_mul(result, n);
        result.data[0] += derivs[k - 1](x)/jet2_factorials[k];
    }
    result = jet2_mul(result, n);
    result.data[0] += f(x);
    return result;
}

static inline jet2
jet2_exp(jet2 a)
{
    static double (*const derivs[2])(double) = {double_exp_deriv1, double_exp_deriv2};
    return jet2_apply(a, double_exp, derivs);
}

static inline jet2
jet2_expm1(jet2 a)
{
    static double (*const derivs[2])(double) = {double_expm1_deriv1, double_expm1_deriv2};
    return jet2_apply(a, double_expm1, derivs);
}

static inline jet2
jet2_log(jet2 a)
{
    static double (*const derivs[2])(double) = {double_log_deriv1, double_log_deriv2};
    return jet2_apply(a, double_log, derivs);
}

static inline jet2
jet2_log1p(jet2 a)
{
    static double (*const derivs[2])(double) = {double_log1p_deriv1, double_log1p_deriv2};
    return jet2_apply(a, double_log1p, derivs);
}

static inline jet2
jet2_sqrt(jet2 a)
{
    static double (*const derivs[2])(double) = {double_sqrt_deriv1, double_sqrt_deriv2};
    return jet2_apply(a, double_sqrt, derivs);
}

static inline jet2
jet2_sin(jet2 a)
{
    static double (*const derivs[2])(double) = {double_sin_deriv1, double_sin_deriv2};
    return jet2_apply(a, double_sin, derivs);
}

static inline jet2
jet2_cos(jet2 a)
{
    static double (*const derivs[2])(double) = {double_cos_deriv1, double_cos_deriv2};
    return jet2_apply(a, double_cos, derivs);
}

static inline jet2
jet2_tan(jet2 a)
{
    static double (*const derivs[2])(double) = {double_tan_deriv1, double_tan_deriv2};
    return jet2_apply(a, double_tan, derivs);
}

static inline jet2
jet2_asin(jet2 a)
{
    static double (*const derivs[2])(double) = {double_asin_deriv1, double_asin_deriv2};
    return jet2_apply(a, double_asin, derivs);
}

static inline jet2
jet2_acos(jet2 a)
{
    static double (*const derivs[2])(double) = {double_acos_deriv1, double_acos_deriv2};
    return jet2_apply(a, double_acos, derivs);
}

static inline jet2
jet2_atan(jet2 a)
{
    static double (*const derivs[2])(double) = {double_atan_deriv1, double_atan_deriv2};
    return jet2_apply(a, double_atan, derivs);
}

static inline jet2
jet2_sinh(jet2 a)
{
    static double (*const derivs[2])(double) = {double_sinh_deriv1, double_sinh_deriv2};
    return jet2_apply(a, double_sinh, derivs);
}

static inline jet2
jet2_cosh(jet2 a)
{
    static double (*const derivs[2])(double) = {double_cosh_deriv1, double_cosh_deriv2};
    return jet2_apply(a, double_cosh, derivs);
}

static inline jet2
jet2_tanh(jet2 a)
{
    static double (*const derivs[2])(double) = {double_tanh_deriv1, double_tanh_deriv2};
    return jet2_apply(a, double_tanh, derivs);
}

static inline jet2
jet2_asinh(jet2 a)
{
    static double (*const derivs[2])(double) = {double_asinh_deriv1, double_asinh_deriv2};
    return jet2_apply(a, double_asinh, derivs);
}

static inline jet2
jet2_acosh(jet2 a)
{
    static double (*const derivs[2])(double) = {double_acosh_deriv1, double_acosh_deriv2};
    return jet2_apply(a, double_acosh, derivs);
}

static inline jet2
jet2_atanh(jet2 a)
{
    static double (*const derivs[2])(double) = {double_atanh_deriv1, double_atanh_deriv2};
    return jet2_apply(a, double_atanh, derivs);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline jet2
jet2_ipow(jet2 a, int64_t b)
{
    jet2 result = jet2_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return jet2_square(a);
    }
    if (b == -1) {
        return jet2_inv(a);
    }
    if (b < 0) {
        a = jet2_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = jet2_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = jet2_square(a);
        }
    }
    return result;
}

/* Real analytic power: requires a positive real part for log(a).
 * Jets are commutative, so lpow and rpow agree. */
static inline jet2
jet2_rpow(jet2 a, jet2 b)
{
    return jet2_exp(jet2_mul(jet2_log(a), b));
}

static inline jet2
jet2_pow(jet2 a, jet2 b)
{
    return jet2_rpow(a, b);
}

static inline jet2
jet2_lpow(jet2 a, jet2 b)
{
    return jet2_rpow(a, b);
}

static inline jet2
jet2_radd(double a, jet2 b)
{
    return jet2_add(jet2_from_re(a), b);
}

static inline jet2
jet2_ladd(jet2 a, double b)
{
    return jet2_add(a, jet2_from_re(b));
}

static inline jet2
jet2_rsub(double a, jet2 b)
{
    return jet2_sub(jet2_from_re(a), b);
}

static inline jet2
jet2_lsub(jet2 a, double b)
{
    return jet2_sub(a, jet2_from_re(b));
}

static inline jet2
jet2_rmul(double a, jet2 b)
{
    return vec3_scale(b, a);
}

static inline jet2
jet2_lmul(jet2 a, double b)
{
    return vec3_scale(a, b);
}

static inline jet2
jet2_rdiv(double a, jet2 b)
{
    return vec3_scale(jet2_inv(b), a);
}

static inline jet2
jet2_ldiv(jet2 a, double b)
{
    return vec3_div_scalar(a, b);
}

static inline jet2
jet2_rpows(double a, jet2 b)
{
    return jet2_rpow(jet2_from_re(a), b);
}

static inline jet2
jet2_lpows(jet2 a, double b)
{
    return jet2_rpow(a, jet2_from_re(b));
}

static inline jet2
jet2_mix(jet2 a, jet2 b, double t)
{
    jet2 delta = jet2_mul(jet2_inv(a), b);
    return jet2_mul(a, jet2_lpows(delta, t));
}

/* |det lmat(a)|^(1/3): the absolute value of the real part. */
static inline double
jet2_norm(jet2 a)
{
    return fabs(a.data[0]);
}

/* Generic minimal polynomial (t - x)^3, ascending coefficients. */
static inline void
jet2_minpoly(jet2 a, double result[4])
{
    double x = a.data[0];

    result[0] = -x*x*x;
    result[1] = 3*x*x;
    result[2] = -3*x;
    result[3] = 1;
}

static inline void
jet2_eivals(jet2 a, complex result[3])
{
    int i;
    for (i = 0; i < 3; i++) {
        result[i] = complex_make((const double[]){a.data[0], 0});
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
jet2_lmat(jet2 a, double result[ALGPY_JET2_N*ALGPY_JET2_N])
{
    jet2 basis = {{0}};
    jet2 product;
    int row, col;
    for (col = 0; col < ALGPY_JET2_N; col++) {
        basis.data[col] = 1;
        product = jet2_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_JET2_N; row++) {
            result[row * ALGPY_JET2_N + col] = product.data[row];
        }
    }
}

static inline void
jet2_rmat(jet2 a, double result[ALGPY_JET2_N*ALGPY_JET2_N])
{
    jet2 basis = {{0}};
    jet2 product;
    int row, col;
    for (row = 0; row < ALGPY_JET2_N; row++) {
        basis.data[row] = 1;
        product = jet2_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_JET2_N; col++) {
            result[row * ALGPY_JET2_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_JET2_H */
