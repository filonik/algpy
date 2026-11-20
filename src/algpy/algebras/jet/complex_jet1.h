/*
 * Complex jets of order 1: truncated polynomials
 * z + z1*e + ... + z1*e^1 in C[e]/(e^2) with complex coefficients,
 * stored as a vec4 of 4 real components: coefficient k is
 * data[2k] + i*data[2k+1]. Over the reals the basis is 1, i, e, i*e, ...,
 * so the generators i and e sit at positions 1 and 2. A holomorphic
 * function applied to z + e gives its normalized complex derivatives
 * f^(k)(z)/k!, k = 0..1, on the principal branches.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX_JET1_H
#define ALGPY_COMPLEX_JET1_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_COMPLEX_JET1_N 4

typedef vec4 complex_jet1;

static inline complex_jet1
complex_jet1_make(const double coeffs[ALGPY_COMPLEX_JET1_N])
{
    complex_jet1 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline complex_jet1
complex_jet1_from_re(double a)
{
    return complex_jet1_make((const double[ALGPY_COMPLEX_JET1_N]){a});
}

static inline double
complex_jet1_to_re(complex_jet1 a)
{
    return a.data[0];
}

/* Complex coefficient k, and its replacement. */
static inline complex
complex_jet1_coefficient(complex_jet1 a, int k)
{
    return complex_make((const double[]){a.data[2*k], a.data[2*k + 1]});
}

static inline void
complex_jet1_set_coefficient(complex_jet1 *a, int k, complex z)
{
    a->data[2*k] = z.data[0];
    a->data[2*k + 1] = z.data[1];
}

static inline complex_jet1
complex_jet1_zero(void)
{
    return vec4_zero();
}

static inline complex_jet1
complex_jet1_one(void)
{
    return complex_jet1_from_re(1);
}

static inline int
complex_jet1_is_real(complex_jet1 a)
{
    return vec4_is_real(a);
}

static inline complex_jet1
complex_jet1_re(complex_jet1 a)
{
    return complex_jet1_from_re(a.data[0]);
}

static inline complex_jet1
complex_jet1_im(complex_jet1 a)
{
    complex_jet1 result = a;
    result.data[0] = 0;
    return result;
}

static inline complex_jet1
complex_jet1_pos(complex_jet1 a)
{
    return vec4_pos(a);
}

static inline complex_jet1
complex_jet1_neg(complex_jet1 a)
{
    return vec4_neg(a);
}

static inline complex_jet1
complex_jet1_add(complex_jet1 a, complex_jet1 b)
{
    return vec4_add(a, b);
}

static inline complex_jet1
complex_jet1_sub(complex_jet1 a, complex_jet1 b)
{
    return vec4_sub(a, b);
}

static inline int
complex_jet1_equal(complex_jet1 a, complex_jet1 b)
{
    return vec4_all_eq(a, b);
}

static inline int
complex_jet1_not_equal(complex_jet1 a, complex_jet1 b)
{
    return !complex_jet1_equal(a, b);
}

/* Truncated Cauchy product of the complex coefficients. */
static inline complex_jet1
complex_jet1_mul(complex_jet1 a, complex_jet1 b)
{
    const double *x = a.data, *y = b.data;
    return complex_jet1_make((const double[]){
        x[0]*y[0] - x[1]*y[1],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] - x[1]*y[3] + x[2]*y[0] - x[3]*y[1],
        x[0]*y[3] + x[1]*y[2] + x[2]*y[1] + x[3]*y[0]});
}

static inline complex_jet1
complex_jet1_square(complex_jet1 a)
{
    return complex_jet1_mul(a, a);
}

/* Complex conjugation i -> -i of every coefficient. */
static inline complex_jet1
complex_jet1_conjugate(complex_jet1 a)
{
    const double *x = a.data;
    return complex_jet1_make((const double[]){x[0], -x[1], x[2], -x[3]});
}

static inline void
complex_jet1_basis(complex_jet1 result[ALGPY_COMPLEX_JET1_N])
{
    int i, j;
    for (i = 0; i < ALGPY_COMPLEX_JET1_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX_JET1_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

/* The generators i and e, at positions 1 and 2. */
static inline void
complex_jet1_gens(complex_jet1 result[2])
{
    result[0] = complex_jet1_make((const double[ALGPY_COMPLEX_JET1_N]){[1] = 1});
    result[1] = complex_jet1_make((const double[ALGPY_COMPLEX_JET1_N]){[2] = 1});
}

static inline void
complex_jet1_coeffs(complex_jet1 result[ALGPY_COMPLEX_JET1_N * ALGPY_COMPLEX_JET1_N])
{
    complex_jet1 basis[ALGPY_COMPLEX_JET1_N];
    int i, j;
    complex_jet1_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX_JET1_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX_JET1_N; j++) {
            result[i * ALGPY_COMPLEX_JET1_N + j] = complex_jet1_mul(basis[i], basis[j]);
        }
    }
}

/* Grade projections: grade k is the complex coefficient of e^k. */
static inline complex_jet1
complex_jet1_grade0(complex_jet1 a)
{
    const double *x = a.data;
    return complex_jet1_make((const double[]){x[0], x[1], 0, 0});
}

static inline complex_jet1
complex_jet1_grade1(complex_jet1 a)
{
    const double *x = a.data;
    return complex_jet1_make((const double[]){0, 0, x[2], x[3]});
}

/* Reciprocal by the series recurrence y0 = 1/z0 and
 * y_k = -(z1*y_(k-1) + ... + zk*y0)/z0. Requires a nonzero leading term. */
static inline complex_jet1
complex_jet1_inv(complex_jet1 a)
{
    complex_jet1 result = {{0}};
    complex r = complex_inv(complex_jet1_coefficient(a, 0));
    int k, i;
    complex_jet1_set_coefficient(&result, 0, r);
    for (k = 1; k < 2; k++) {
        complex sum = complex_zero();
        for (i = 1; i <= k; i++) {
            sum = complex_add(sum, complex_mul(complex_jet1_coefficient(a, i),
                                               complex_jet1_coefficient(result, k - i)));
        }
        complex_jet1_set_coefficient(&result, k, complex_neg(complex_mul(sum, r)));
    }
    return result;
}

static inline complex_jet1
complex_jet1_div(complex_jet1 a, complex_jet1 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return complex_jet1_mul(a, complex_jet1_inv(b));
}

static const double complex_jet1_factorials[2] = {1., 1.};

/*
 * Taylor expansion about the leading coefficient z: with a = z + n,
 * n nilpotent, f(a) = sum over k of f^(k)(z)/k! * n^k, by Horner's rule.
 * derivs[k-1] is the k-th complex derivative on the principal branch.
 */
static inline complex_jet1
complex_jet1_apply(complex_jet1 a, complex (*f)(complex), complex (*const derivs[1])(complex))
{
    complex z = complex_jet1_coefficient(a, 0), c;
    complex_jet1 n = a, result = {{0}};
    int k;
    n.data[0] = n.data[1] = 0;
    complex_jet1_set_coefficient(&result, 0, vec2_div_scalar(derivs[0](z), complex_jet1_factorials[1]));
    for (k = 0; k >= 1; k--) {
        result = complex_jet1_mul(result, n);
        c = vec2_div_scalar(derivs[k - 1](z), complex_jet1_factorials[k]);
        result.data[0] += c.data[0];
        result.data[1] += c.data[1];
    }
    result = complex_jet1_mul(result, n);
    c = f(z);
    result.data[0] += c.data[0];
    result.data[1] += c.data[1];
    return result;
}

static inline complex_jet1
complex_jet1_exp(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_exp_deriv1};
    return complex_jet1_apply(a, complex_exp, derivs);
}

static inline complex_jet1
complex_jet1_expm1(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_expm1_deriv1};
    return complex_jet1_apply(a, complex_expm1, derivs);
}

static inline complex_jet1
complex_jet1_log(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_log_deriv1};
    return complex_jet1_apply(a, complex_log, derivs);
}

static inline complex_jet1
complex_jet1_log1p(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_log1p_deriv1};
    return complex_jet1_apply(a, complex_log1p, derivs);
}

static inline complex_jet1
complex_jet1_sqrt(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_sqrt_deriv1};
    return complex_jet1_apply(a, complex_sqrt, derivs);
}

static inline complex_jet1
complex_jet1_sin(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_sin_deriv1};
    return complex_jet1_apply(a, complex_sin, derivs);
}

static inline complex_jet1
complex_jet1_cos(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_cos_deriv1};
    return complex_jet1_apply(a, complex_cos, derivs);
}

static inline complex_jet1
complex_jet1_tan(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_tan_deriv1};
    return complex_jet1_apply(a, complex_tan, derivs);
}

static inline complex_jet1
complex_jet1_asin(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_asin_deriv1};
    return complex_jet1_apply(a, complex_asin, derivs);
}

static inline complex_jet1
complex_jet1_acos(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_acos_deriv1};
    return complex_jet1_apply(a, complex_acos, derivs);
}

static inline complex_jet1
complex_jet1_atan(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_atan_deriv1};
    return complex_jet1_apply(a, complex_atan, derivs);
}

static inline complex_jet1
complex_jet1_sinh(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_sinh_deriv1};
    return complex_jet1_apply(a, complex_sinh, derivs);
}

static inline complex_jet1
complex_jet1_cosh(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_cosh_deriv1};
    return complex_jet1_apply(a, complex_cosh, derivs);
}

static inline complex_jet1
complex_jet1_tanh(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_tanh_deriv1};
    return complex_jet1_apply(a, complex_tanh, derivs);
}

static inline complex_jet1
complex_jet1_asinh(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_asinh_deriv1};
    return complex_jet1_apply(a, complex_asinh, derivs);
}

static inline complex_jet1
complex_jet1_acosh(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_acosh_deriv1};
    return complex_jet1_apply(a, complex_acosh, derivs);
}

static inline complex_jet1
complex_jet1_atanh(complex_jet1 a)
{
    static complex (*const derivs[1])(complex) = {complex_atanh_deriv1};
    return complex_jet1_apply(a, complex_atanh, derivs);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline complex_jet1
complex_jet1_ipow(complex_jet1 a, int64_t b)
{
    complex_jet1 result = complex_jet1_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return complex_jet1_square(a);
    }
    if (b == -1) {
        return complex_jet1_inv(a);
    }
    if (b < 0) {
        a = complex_jet1_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = complex_jet1_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = complex_jet1_square(a);
        }
    }
    return result;
}

/* Principal power exp(log(a)*b); complex jets commute, so lpow = rpow. */
static inline complex_jet1
complex_jet1_rpow(complex_jet1 a, complex_jet1 b)
{
    return complex_jet1_exp(complex_jet1_mul(complex_jet1_log(a), b));
}

static inline complex_jet1
complex_jet1_pow(complex_jet1 a, complex_jet1 b)
{
    return complex_jet1_rpow(a, b);
}

static inline complex_jet1
complex_jet1_lpow(complex_jet1 a, complex_jet1 b)
{
    return complex_jet1_rpow(a, b);
}

static inline complex_jet1
complex_jet1_radd(double a, complex_jet1 b)
{
    return complex_jet1_add(complex_jet1_from_re(a), b);
}

static inline complex_jet1
complex_jet1_ladd(complex_jet1 a, double b)
{
    return complex_jet1_add(a, complex_jet1_from_re(b));
}

static inline complex_jet1
complex_jet1_rsub(double a, complex_jet1 b)
{
    return complex_jet1_sub(complex_jet1_from_re(a), b);
}

static inline complex_jet1
complex_jet1_lsub(complex_jet1 a, double b)
{
    return complex_jet1_sub(a, complex_jet1_from_re(b));
}

static inline complex_jet1
complex_jet1_rmul(double a, complex_jet1 b)
{
    return vec4_scale(b, a);
}

static inline complex_jet1
complex_jet1_lmul(complex_jet1 a, double b)
{
    return vec4_scale(a, b);
}

static inline complex_jet1
complex_jet1_rdiv(double a, complex_jet1 b)
{
    return vec4_scale(complex_jet1_inv(b), a);
}

static inline complex_jet1
complex_jet1_ldiv(complex_jet1 a, double b)
{
    return vec4_div_scalar(a, b);
}

static inline complex_jet1
complex_jet1_rpows(double a, complex_jet1 b)
{
    return complex_jet1_rpow(complex_jet1_from_re(a), b);
}

static inline complex_jet1
complex_jet1_lpows(complex_jet1 a, double b)
{
    return complex_jet1_rpow(a, complex_jet1_from_re(b));
}

static inline complex_jet1
complex_jet1_mix(complex_jet1 a, complex_jet1 b, double t)
{
    complex_jet1 delta = complex_jet1_mul(complex_jet1_inv(a), b);
    return complex_jet1_mul(a, complex_jet1_lpows(delta, t));
}

/* |det lmat(a)|^(1/4): the modulus of the leading coefficient. */
static inline double
complex_jet1_norm(complex_jet1 a)
{
    return hypot(a.data[0], a.data[1]);
}

/* Generic minimal polynomial over the reals, ascending:
 * ((t - z)(t - conj(z)))^2 = (t^2 - 2 Re(z) t + |z|^2)^2. */
static inline void
complex_jet1_minpoly(complex_jet1 a, double result[5])
{
    double q[3] = {a.data[0]*a.data[0] + a.data[1]*a.data[1], -2*a.data[0], 1};
    int i, j, degree = 0;
    for (i = 0; i <= 4; i++) {
        result[i] = i == 0;
    }
    for (i = 0; i < 2; i++) {
        for (j = degree + 2; j >= 0; j--) {
            double value = q[0]*result[j];
            if (j >= 1) {
                value += q[1]*result[j - 1];
            }
            if (j >= 2) {
                value += q[2]*result[j - 2];
            }
            result[j] = value;
        }
        degree += 2;
    }
}

/* The leading coefficient and its conjugate, each 2 times. */
static inline void
complex_jet1_eivals(complex_jet1 a, complex result[4])
{
    int i;
    for (i = 0; i < 2; i++) {
        result[2*i] = complex_make((const double[]){a.data[0], a.data[1]});
        result[2*i + 1] = complex_make((const double[]){a.data[0], -a.data[1]});
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
complex_jet1_lmat(complex_jet1 a, double result[ALGPY_COMPLEX_JET1_N*ALGPY_COMPLEX_JET1_N])
{
    complex_jet1 basis = {{0}};
    complex_jet1 product;
    int row, col;
    for (col = 0; col < ALGPY_COMPLEX_JET1_N; col++) {
        basis.data[col] = 1;
        product = complex_jet1_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_COMPLEX_JET1_N; row++) {
            result[row * ALGPY_COMPLEX_JET1_N + col] = product.data[row];
        }
    }
}

static inline void
complex_jet1_rmat(complex_jet1 a, double result[ALGPY_COMPLEX_JET1_N*ALGPY_COMPLEX_JET1_N])
{
    complex_jet1 basis = {{0}};
    complex_jet1 product;
    int row, col;
    for (row = 0; row < ALGPY_COMPLEX_JET1_N; row++) {
        basis.data[row] = 1;
        product = complex_jet1_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_COMPLEX_JET1_N; col++) {
            result[row * ALGPY_COMPLEX_JET1_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_COMPLEX_JET1_H */
