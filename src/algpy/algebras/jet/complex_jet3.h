/*
 * Complex jets of order 3: truncated polynomials
 * z + z1*e + ... + z3*e^3 in C[e]/(e^4) with complex coefficients,
 * stored as a vec8 of 8 real components: coefficient k is
 * data[2k] + i*data[2k+1]. Over the reals the basis is 1, i, e, i*e, ...,
 * so the generators i and e sit at positions 1 and 2. A holomorphic
 * function applied to z + e gives its normalized complex derivatives
 * f^(k)(z)/k!, k = 0..3, on the principal branches.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX_JET3_H
#define ALGPY_COMPLEX_JET3_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../vectorspaces/vec8.h"

#define ALGPY_COMPLEX_JET3_N 8

typedef vec8 complex_jet3;

static inline complex_jet3
complex_jet3_make(const double coeffs[ALGPY_COMPLEX_JET3_N])
{
    complex_jet3 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline complex_jet3
complex_jet3_from_re(double a)
{
    return complex_jet3_make((const double[ALGPY_COMPLEX_JET3_N]){a});
}

static inline double
complex_jet3_to_re(complex_jet3 a)
{
    return a.data[0];
}

/* Complex coefficient k, and its replacement. */
static inline complex
complex_jet3_coefficient(complex_jet3 a, int k)
{
    return complex_make((const double[]){a.data[2*k], a.data[2*k + 1]});
}

static inline void
complex_jet3_set_coefficient(complex_jet3 *a, int k, complex z)
{
    a->data[2*k] = z.data[0];
    a->data[2*k + 1] = z.data[1];
}

static inline complex_jet3
complex_jet3_zero(void)
{
    return vec8_zero();
}

static inline complex_jet3
complex_jet3_one(void)
{
    return complex_jet3_from_re(1);
}

static inline int
complex_jet3_is_real(complex_jet3 a)
{
    return vec8_is_real(a);
}

static inline complex_jet3
complex_jet3_re(complex_jet3 a)
{
    return complex_jet3_from_re(a.data[0]);
}

static inline complex_jet3
complex_jet3_im(complex_jet3 a)
{
    complex_jet3 result = a;
    result.data[0] = 0;
    return result;
}

static inline complex_jet3
complex_jet3_pos(complex_jet3 a)
{
    return vec8_pos(a);
}

static inline complex_jet3
complex_jet3_neg(complex_jet3 a)
{
    return vec8_neg(a);
}

static inline complex_jet3
complex_jet3_add(complex_jet3 a, complex_jet3 b)
{
    return vec8_add(a, b);
}

static inline complex_jet3
complex_jet3_sub(complex_jet3 a, complex_jet3 b)
{
    return vec8_sub(a, b);
}

static inline int
complex_jet3_equal(complex_jet3 a, complex_jet3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
complex_jet3_not_equal(complex_jet3 a, complex_jet3 b)
{
    return !complex_jet3_equal(a, b);
}

/* Truncated Cauchy product of the complex coefficients. */
static inline complex_jet3
complex_jet3_mul(complex_jet3 a, complex_jet3 b)
{
    const double *x = a.data, *y = b.data;
    return complex_jet3_make((const double[]){
        x[0]*y[0] - x[1]*y[1],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] - x[1]*y[3] + x[2]*y[0] - x[3]*y[1],
        x[0]*y[3] + x[1]*y[2] + x[2]*y[1] + x[3]*y[0],
        x[0]*y[4] - x[1]*y[5] + x[2]*y[2] - x[3]*y[3] + x[4]*y[0] - x[5]*y[1],
        x[0]*y[5] + x[1]*y[4] + x[2]*y[3] + x[3]*y[2] + x[4]*y[1] + x[5]*y[0],
        x[0]*y[6] - x[1]*y[7] + x[2]*y[4] - x[3]*y[5] + x[4]*y[2] - x[5]*y[3] + x[6]*y[0] - x[7]*y[1],
        x[0]*y[7] + x[1]*y[6] + x[2]*y[5] + x[3]*y[4] + x[4]*y[3] + x[5]*y[2] + x[6]*y[1] + x[7]*y[0]});
}

static inline complex_jet3
complex_jet3_square(complex_jet3 a)
{
    return complex_jet3_mul(a, a);
}

/* Complex conjugation i -> -i of every coefficient. */
static inline complex_jet3
complex_jet3_conjugate(complex_jet3 a)
{
    const double *x = a.data;
    return complex_jet3_make((const double[]){x[0], -x[1], x[2], -x[3], x[4], -x[5], x[6], -x[7]});
}

static inline void
complex_jet3_basis(complex_jet3 result[ALGPY_COMPLEX_JET3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_COMPLEX_JET3_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX_JET3_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

/* The generators i and e, at positions 1 and 2. */
static inline void
complex_jet3_gens(complex_jet3 result[2])
{
    result[0] = complex_jet3_make((const double[ALGPY_COMPLEX_JET3_N]){[1] = 1});
    result[1] = complex_jet3_make((const double[ALGPY_COMPLEX_JET3_N]){[2] = 1});
}

static inline void
complex_jet3_coeffs(complex_jet3 result[ALGPY_COMPLEX_JET3_N * ALGPY_COMPLEX_JET3_N])
{
    complex_jet3 basis[ALGPY_COMPLEX_JET3_N];
    int i, j;
    complex_jet3_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX_JET3_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX_JET3_N; j++) {
            result[i * ALGPY_COMPLEX_JET3_N + j] = complex_jet3_mul(basis[i], basis[j]);
        }
    }
}

/* Grade projections: grade k is the complex coefficient of e^k. */
static inline complex_jet3
complex_jet3_grade0(complex_jet3 a)
{
    const double *x = a.data;
    return complex_jet3_make((const double[]){x[0], x[1], 0, 0, 0, 0, 0, 0});
}

static inline complex_jet3
complex_jet3_grade1(complex_jet3 a)
{
    const double *x = a.data;
    return complex_jet3_make((const double[]){0, 0, x[2], x[3], 0, 0, 0, 0});
}

static inline complex_jet3
complex_jet3_grade2(complex_jet3 a)
{
    const double *x = a.data;
    return complex_jet3_make((const double[]){0, 0, 0, 0, x[4], x[5], 0, 0});
}

static inline complex_jet3
complex_jet3_grade3(complex_jet3 a)
{
    const double *x = a.data;
    return complex_jet3_make((const double[]){0, 0, 0, 0, 0, 0, x[6], x[7]});
}

/* Reciprocal by the series recurrence y0 = 1/z0 and
 * y_k = -(z1*y_(k-1) + ... + zk*y0)/z0. Requires a nonzero leading term. */
static inline complex_jet3
complex_jet3_inv(complex_jet3 a)
{
    complex_jet3 result = {{0}};
    complex r = complex_inv(complex_jet3_coefficient(a, 0));
    int k, i;
    complex_jet3_set_coefficient(&result, 0, r);
    for (k = 1; k < 4; k++) {
        complex sum = complex_zero();
        for (i = 1; i <= k; i++) {
            sum = complex_add(sum, complex_mul(complex_jet3_coefficient(a, i),
                                               complex_jet3_coefficient(result, k - i)));
        }
        complex_jet3_set_coefficient(&result, k, complex_neg(complex_mul(sum, r)));
    }
    return result;
}

static inline complex_jet3
complex_jet3_div(complex_jet3 a, complex_jet3 b)
{
    if (vec8_is_real(b)) {
        return vec8_div_scalar(a, b.data[0]);
    }
    return complex_jet3_mul(a, complex_jet3_inv(b));
}

static const double complex_jet3_factorials[4] = {1., 1., 2., 6.};

/*
 * Taylor expansion about the leading coefficient z: with a = z + n,
 * n nilpotent, f(a) = sum over k of f^(k)(z)/k! * n^k, by Horner's rule.
 * derivs[k-1] is the k-th complex derivative on the principal branch.
 */
static inline complex_jet3
complex_jet3_apply(complex_jet3 a, complex (*f)(complex), complex (*const derivs[3])(complex))
{
    complex z = complex_jet3_coefficient(a, 0), c;
    complex_jet3 n = a, result = {{0}};
    int k;
    n.data[0] = n.data[1] = 0;
    complex_jet3_set_coefficient(&result, 0, vec2_div_scalar(derivs[2](z), complex_jet3_factorials[3]));
    for (k = 2; k >= 1; k--) {
        result = complex_jet3_mul(result, n);
        c = vec2_div_scalar(derivs[k - 1](z), complex_jet3_factorials[k]);
        result.data[0] += c.data[0];
        result.data[1] += c.data[1];
    }
    result = complex_jet3_mul(result, n);
    c = f(z);
    result.data[0] += c.data[0];
    result.data[1] += c.data[1];
    return result;
}

static inline complex_jet3
complex_jet3_exp(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_exp_deriv1, complex_exp_deriv2, complex_exp_deriv3};
    return complex_jet3_apply(a, complex_exp, derivs);
}

static inline complex_jet3
complex_jet3_expm1(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_expm1_deriv1, complex_expm1_deriv2, complex_expm1_deriv3};
    return complex_jet3_apply(a, complex_expm1, derivs);
}

static inline complex_jet3
complex_jet3_log(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_log_deriv1, complex_log_deriv2, complex_log_deriv3};
    return complex_jet3_apply(a, complex_log, derivs);
}

static inline complex_jet3
complex_jet3_log1p(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_log1p_deriv1, complex_log1p_deriv2, complex_log1p_deriv3};
    return complex_jet3_apply(a, complex_log1p, derivs);
}

static inline complex_jet3
complex_jet3_sqrt(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_sqrt_deriv1, complex_sqrt_deriv2, complex_sqrt_deriv3};
    return complex_jet3_apply(a, complex_sqrt, derivs);
}

static inline complex_jet3
complex_jet3_sin(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_sin_deriv1, complex_sin_deriv2, complex_sin_deriv3};
    return complex_jet3_apply(a, complex_sin, derivs);
}

static inline complex_jet3
complex_jet3_cos(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_cos_deriv1, complex_cos_deriv2, complex_cos_deriv3};
    return complex_jet3_apply(a, complex_cos, derivs);
}

static inline complex_jet3
complex_jet3_tan(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_tan_deriv1, complex_tan_deriv2, complex_tan_deriv3};
    return complex_jet3_apply(a, complex_tan, derivs);
}

static inline complex_jet3
complex_jet3_asin(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_asin_deriv1, complex_asin_deriv2, complex_asin_deriv3};
    return complex_jet3_apply(a, complex_asin, derivs);
}

static inline complex_jet3
complex_jet3_acos(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_acos_deriv1, complex_acos_deriv2, complex_acos_deriv3};
    return complex_jet3_apply(a, complex_acos, derivs);
}

static inline complex_jet3
complex_jet3_atan(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_atan_deriv1, complex_atan_deriv2, complex_atan_deriv3};
    return complex_jet3_apply(a, complex_atan, derivs);
}

static inline complex_jet3
complex_jet3_sinh(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_sinh_deriv1, complex_sinh_deriv2, complex_sinh_deriv3};
    return complex_jet3_apply(a, complex_sinh, derivs);
}

static inline complex_jet3
complex_jet3_cosh(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_cosh_deriv1, complex_cosh_deriv2, complex_cosh_deriv3};
    return complex_jet3_apply(a, complex_cosh, derivs);
}

static inline complex_jet3
complex_jet3_tanh(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_tanh_deriv1, complex_tanh_deriv2, complex_tanh_deriv3};
    return complex_jet3_apply(a, complex_tanh, derivs);
}

static inline complex_jet3
complex_jet3_asinh(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_asinh_deriv1, complex_asinh_deriv2, complex_asinh_deriv3};
    return complex_jet3_apply(a, complex_asinh, derivs);
}

static inline complex_jet3
complex_jet3_acosh(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_acosh_deriv1, complex_acosh_deriv2, complex_acosh_deriv3};
    return complex_jet3_apply(a, complex_acosh, derivs);
}

static inline complex_jet3
complex_jet3_atanh(complex_jet3 a)
{
    static complex (*const derivs[3])(complex) = {complex_atanh_deriv1, complex_atanh_deriv2, complex_atanh_deriv3};
    return complex_jet3_apply(a, complex_atanh, derivs);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline complex_jet3
complex_jet3_ipow(complex_jet3 a, int64_t b)
{
    complex_jet3 result = complex_jet3_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return complex_jet3_square(a);
    }
    if (b == -1) {
        return complex_jet3_inv(a);
    }
    if (b < 0) {
        a = complex_jet3_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = complex_jet3_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = complex_jet3_square(a);
        }
    }
    return result;
}

/* Principal power exp(log(a)*b); complex jets commute, so lpow = rpow. */
static inline complex_jet3
complex_jet3_rpow(complex_jet3 a, complex_jet3 b)
{
    return complex_jet3_exp(complex_jet3_mul(complex_jet3_log(a), b));
}

static inline complex_jet3
complex_jet3_pow(complex_jet3 a, complex_jet3 b)
{
    return complex_jet3_rpow(a, b);
}

static inline complex_jet3
complex_jet3_lpow(complex_jet3 a, complex_jet3 b)
{
    return complex_jet3_rpow(a, b);
}

static inline complex_jet3
complex_jet3_radd(double a, complex_jet3 b)
{
    return complex_jet3_add(complex_jet3_from_re(a), b);
}

static inline complex_jet3
complex_jet3_ladd(complex_jet3 a, double b)
{
    return complex_jet3_add(a, complex_jet3_from_re(b));
}

static inline complex_jet3
complex_jet3_rsub(double a, complex_jet3 b)
{
    return complex_jet3_sub(complex_jet3_from_re(a), b);
}

static inline complex_jet3
complex_jet3_lsub(complex_jet3 a, double b)
{
    return complex_jet3_sub(a, complex_jet3_from_re(b));
}

static inline complex_jet3
complex_jet3_rmul(double a, complex_jet3 b)
{
    return vec8_scale(b, a);
}

static inline complex_jet3
complex_jet3_lmul(complex_jet3 a, double b)
{
    return vec8_scale(a, b);
}

static inline complex_jet3
complex_jet3_rdiv(double a, complex_jet3 b)
{
    return vec8_scale(complex_jet3_inv(b), a);
}

static inline complex_jet3
complex_jet3_ldiv(complex_jet3 a, double b)
{
    return vec8_div_scalar(a, b);
}

static inline complex_jet3
complex_jet3_rpows(double a, complex_jet3 b)
{
    return complex_jet3_rpow(complex_jet3_from_re(a), b);
}

static inline complex_jet3
complex_jet3_lpows(complex_jet3 a, double b)
{
    return complex_jet3_rpow(a, complex_jet3_from_re(b));
}

static inline complex_jet3
complex_jet3_mix(complex_jet3 a, complex_jet3 b, double t)
{
    complex_jet3 delta = complex_jet3_mul(complex_jet3_inv(a), b);
    return complex_jet3_mul(a, complex_jet3_lpows(delta, t));
}

/* |det lmat(a)|^(1/8): the modulus of the leading coefficient. */
static inline double
complex_jet3_norm(complex_jet3 a)
{
    return hypot(a.data[0], a.data[1]);
}

/* Generic minimal polynomial over the reals, ascending:
 * ((t - z)(t - conj(z)))^4 = (t^2 - 2 Re(z) t + |z|^2)^4. */
static inline void
complex_jet3_minpoly(complex_jet3 a, double result[9])
{
    double q[3] = {a.data[0]*a.data[0] + a.data[1]*a.data[1], -2*a.data[0], 1};
    int i, j, degree = 0;
    for (i = 0; i <= 8; i++) {
        result[i] = i == 0;
    }
    for (i = 0; i < 4; i++) {
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

/* The leading coefficient and its conjugate, each 4 times. */
static inline void
complex_jet3_eivals(complex_jet3 a, complex result[8])
{
    int i;
    for (i = 0; i < 4; i++) {
        result[2*i] = complex_make((const double[]){a.data[0], a.data[1]});
        result[2*i + 1] = complex_make((const double[]){a.data[0], -a.data[1]});
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
complex_jet3_lmat(complex_jet3 a, double result[ALGPY_COMPLEX_JET3_N*ALGPY_COMPLEX_JET3_N])
{
    complex_jet3 basis = {{0}};
    complex_jet3 product;
    int row, col;
    for (col = 0; col < ALGPY_COMPLEX_JET3_N; col++) {
        basis.data[col] = 1;
        product = complex_jet3_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_COMPLEX_JET3_N; row++) {
            result[row * ALGPY_COMPLEX_JET3_N + col] = product.data[row];
        }
    }
}

static inline void
complex_jet3_rmat(complex_jet3 a, double result[ALGPY_COMPLEX_JET3_N*ALGPY_COMPLEX_JET3_N])
{
    complex_jet3 basis = {{0}};
    complex_jet3 product;
    int row, col;
    for (row = 0; row < ALGPY_COMPLEX_JET3_N; row++) {
        basis.data[row] = 1;
        product = complex_jet3_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_COMPLEX_JET3_N; col++) {
            result[row * ALGPY_COMPLEX_JET3_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_COMPLEX_JET3_H */
