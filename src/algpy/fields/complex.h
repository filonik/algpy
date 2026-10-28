/*
 * Shared arithmetic of complex numbers (x + y*i), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX_H
#define ALGPY_COMPLEX_H

#include <math.h>
#include <complex.h>
#include <string.h>

#include "../vectorspaces/vec2.h"

#define ALGPY_COMPLEX_N ALGPY_VEC2_N

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* C99 <complex.h> defines complex as a macro; reserve this name for vec2. */
#undef complex

typedef vec2 complex;

static inline complex
complex_make(double x, double y)
{
    complex r = {{x, y}};
    return r;
}

static inline complex
complex_re(complex a)
{
    return complex_make(a.data[0], 0);
}

static inline complex
complex_im(complex a)
{
    return complex_make(0, a.data[1]);
}

static inline complex
complex_from_re(double a)
{
    return complex_make(a, 0);
}

static inline double
complex_to_re(complex a)
{
    return a.data[0];
}

static inline complex
complex_one(void)
{
    return complex_from_re(1);
}

static inline complex
complex_mul(complex x, complex y)
{
    return complex_make(
        x.data[0]*y.data[0] - x.data[1]*y.data[1],
        x.data[0]*y.data[1] + x.data[1]*y.data[0]);
}

static inline void
complex_basis(complex result[ALGPY_COMPLEX_N])
{
    int i, j;
    for (i = 0; i < ALGPY_COMPLEX_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
complex_coeffs(
    complex result[ALGPY_COMPLEX_N * ALGPY_COMPLEX_N])
{
    complex basis[ALGPY_COMPLEX_N];
    int i, j;
    complex_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX_N; j++) {
            result[i * ALGPY_COMPLEX_N + j] =
                complex_mul(basis[i], basis[j]);
        }
    }
}

static inline complex
complex_conjugate(complex a)
{
    return complex_make(a.data[0], -a.data[1]);
}

/* Norm a*conj(a) = x^2 + y^2 */
static inline double
complex_norm2(complex a)
{
    return a.data[0]*a.data[0] + a.data[1]*a.data[1];
}

static inline complex
complex_inv(complex a)
{
    double n = complex_norm2(a);
    return complex_make(a.data[0] / n, -a.data[1] / n);
}

static inline complex
complex_div(complex a, complex b)
{
    return complex_mul(a, complex_inv(b));
}

static inline complex
complex_exp(complex a)
{
    double e = exp(a.data[0]);
    return complex_make(e * cos(a.data[1]), e * sin(a.data[1]));
}

static inline complex
complex_log(complex a)
{
    double r = hypot(a.data[0], a.data[1]);
    return complex_make(r == 0 ? -INFINITY : log(r),
                         atan2(a.data[1], a.data[0]));
}

static inline complex
complex_pow(complex a, complex b)
{
    return complex_exp(complex_mul(complex_log(a), b));
}

static inline complex
complex_radd(double a, complex b) {
  return vec2_add(complex_from_re(a), b);
}

static inline complex
complex_ladd(complex a, double b) {
  return vec2_add(a, complex_from_re(b));
}

static inline complex
complex_rsub(double a, complex b) {
  return vec2_sub(complex_from_re(a), b);
}

static inline complex
complex_lsub(complex a, double b) {
  return vec2_sub(a, complex_from_re(b));
}

static inline complex
complex_rmul(double a, complex b) {
  return complex_mul(complex_from_re(a), b);
}

static inline complex
complex_lmul(complex a, double b) {
  return complex_mul(a, complex_from_re(b));
}

static inline complex
complex_rdiv(double a, complex b) {
  return complex_div(complex_from_re(a), b);
}

static inline complex
complex_ldiv(complex a, double b) {
  return complex_div(a, complex_from_re(b));
}

static inline complex
complex_rpow(double a, complex b) {
  return complex_pow(complex_from_re(a), b);
}

static inline complex
complex_lpow(complex a, double b) {
  return complex_pow(a, complex_from_re(b));
}

static inline complex
complex_mix(complex a, complex b, double t)
{
    complex delta = complex_mul(complex_inv(a), b);
    return complex_mul(a, complex_lpow(delta, t));
}

static inline complex
complex_expm1(complex a)
{
    return complex_lsub(complex_exp(a), 1);
}

static inline complex
complex_log1p(complex a)
{
    return complex_log(complex_ladd(a, 1));
}

static inline complex
complex_sqrt(complex a)
{
    double r = hypot(a.data[0], a.data[1]);
    double re = sqrt(fmax(0, (r + a.data[0]) / 2));
    double im = sqrt(fmax(0, (r - a.data[0]) / 2));
    return complex_make(re, copysign(im, a.data[1]));
}

static inline complex
complex_cos(complex a)
{
    return complex_make(cos(a.data[0]) * cosh(a.data[1]),
                         -sin(a.data[0]) * sinh(a.data[1]));
}

static inline complex
complex_sin(complex a)
{
    return complex_make(sin(a.data[0]) * cosh(a.data[1]),
                         cos(a.data[0]) * sinh(a.data[1]));
}

static inline complex
complex_tan(complex a)
{
    return complex_div(complex_sin(a), complex_cos(a));
}

static inline complex
complex_acos(complex a)
{
    double _Complex z;
    /* C99 complex has the layout of a double[2]. Copying preserves signed
     * zeros on branch cuts without arithmetic on the components. */
    memcpy(&z, a.data, sizeof(z));
    z = cacos(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_asin(complex a)
{
    double _Complex z;
    memcpy(&z, a.data, sizeof(z));
    z = casin(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_atan(complex a)
{
    double _Complex z;
    memcpy(&z, a.data, sizeof(z));
    z = catan(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_cosh(complex a)
{
    return complex_make(cosh(a.data[0]) * cos(a.data[1]),
                         sinh(a.data[0]) * sin(a.data[1]));
}

static inline complex
complex_sinh(complex a)
{
    return complex_make(sinh(a.data[0]) * cos(a.data[1]),
                         cosh(a.data[0]) * sin(a.data[1]));
}

static inline complex
complex_tanh(complex a)
{
    return complex_div(complex_sinh(a), complex_cosh(a));
}

static inline complex
complex_acosh(complex a)
{
    double _Complex z;
    memcpy(&z, a.data, sizeof(z));
    z = cacosh(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_asinh(complex a)
{
    double _Complex z;
    memcpy(&z, a.data, sizeof(z));
    z = casinh(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_atanh(complex a)
{
    double _Complex z;
    memcpy(&z, a.data, sizeof(z));
    z = catanh(z);
    return complex_make(creal(z), cimag(z));
}

static inline void
complex_minpoly(complex a, double result[3])
{
    result[0] = complex_norm2(a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

static inline void
complex_eivals(complex a, complex result[2])
{
    result[0] = complex_conjugate(a);
    result[1] = a;
}

/* Square root of a real number: sqrt(k) or i*sqrt(-k). */
static inline complex
complex_sqrt_real(double k)
{
    if (k < 0) {
        return complex_make(0, sqrt(-k));
    }
    return complex_make(sqrt(k), 0);
}

/* Principal d-th root of a real number: |k|^(1/d) times exp(i*pi/d) if k < 0. */
static inline complex
complex_root_real(double k, int d)
{
    if (k < 0 && d == 2) {
        return complex_make(0, sqrt(-k));
    }
    if (k < 0) {
        double r = pow(-k, 1. / d);
        double t = M_PI / d;
        return complex_make(r * cos(t), r * sin(t));
    }
    return complex_make(pow(k, 1. / d), 0);
}

/* Common algebra API. */
static inline complex
complex_zero(void)
{
    return vec2_zero();
}

static inline int
complex_is_real(complex a)
{
    return vec2_is_real(a);
}

static inline complex
complex_pos(complex a)
{
    return vec2_pos(a);
}

static inline complex
complex_neg(complex a)
{
    return vec2_neg(a);
}

static inline complex
complex_add(complex a, complex b)
{
    return vec2_add(a, b);
}

static inline complex
complex_sub(complex a, complex b)
{
    return vec2_sub(a, b);
}

static inline int
complex_equal(complex a, complex b)
{
    return vec2_all_eq(a, b);
}

static inline int
complex_not_equal(complex a, complex b)
{
    return !complex_equal(a, b);
}

static inline complex
complex_square(complex a)
{
    return complex_mul(a, a);
}

static inline double
complex_norm(complex a)
{
    return sqrt(complex_norm2(a));
}

/* Ordinary complex derivatives f^(N), not Taylor coefficients divided by N!.
 * Branches agree with the principal field functions. Values on a branch cut
 * are the corresponding one-sided limits, selected by signed zeros.
 * Singular points follow IEEE floating-point propagation. */
static inline complex
complex_exp_deriv1(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_exp_deriv2(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_exp_deriv3(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_exp_deriv4(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_expm1_deriv1(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_expm1_deriv2(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_expm1_deriv3(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_expm1_deriv4(complex a)
{
    return complex_exp(a);
}

static inline complex
complex_log_deriv1(complex a)
{
    return complex_inv(a);
}

static inline complex
complex_log_deriv2(complex a)
{
    return complex_neg(complex_square(complex_inv(a)));
}

static inline complex
complex_log_deriv3(complex a)
{
    complex r = complex_inv(a);
    return vec2_scale(complex_mul(complex_square(r), r), 2);
}

static inline complex
complex_log_deriv4(complex a)
{
    complex r2 = complex_square(complex_inv(a));
    return vec2_scale(complex_square(r2), -6);
}

static inline complex
complex_log1p_deriv1(complex a)
{
    return complex_log_deriv1(complex_make(1+a.data[0], a.data[1]));
}

static inline complex
complex_log1p_deriv2(complex a)
{
    return complex_log_deriv2(complex_make(1+a.data[0], a.data[1]));
}

static inline complex
complex_log1p_deriv3(complex a)
{
    return complex_log_deriv3(complex_make(1+a.data[0], a.data[1]));
}

static inline complex
complex_log1p_deriv4(complex a)
{
    return complex_log_deriv4(complex_make(1+a.data[0], a.data[1]));
}

static inline complex
complex_inv_deriv1(complex a)
{
    return complex_neg(complex_square(complex_inv(a)));
}

static inline complex
complex_inv_deriv2(complex a)
{
    complex r = complex_inv(a);
    return vec2_scale(complex_mul(complex_square(r), r), 2);
}

static inline complex
complex_inv_deriv3(complex a)
{
    complex r2 = complex_square(complex_inv(a));
    return vec2_scale(complex_square(r2), -6);
}

static inline complex
complex_inv_deriv4(complex a)
{
    complex r = complex_inv(a);
    return vec2_scale(complex_mul(complex_square(complex_square(r)), r), 24);
}

static inline complex
complex_square_deriv1(complex a)
{
    return vec2_scale(a, 2);
}

static inline complex
complex_square_deriv2(complex a)
{
    (void)a;
    return complex_from_re(2);
}

static inline complex
complex_square_deriv3(complex a)
{
    (void)a;
    return complex_zero();
}

static inline complex
complex_square_deriv4(complex a)
{
    (void)a;
    return complex_zero();
}

static inline complex
complex_sqrt_deriv1(complex a)
{
    return vec2_scale(complex_inv(complex_sqrt(a)), 0.5);
}

static inline complex
complex_sqrt_deriv2(complex a)
{
    return vec2_scale(complex_mul(complex_inv(a), complex_inv(complex_sqrt(a))), -0.25);
}

static inline complex
complex_sqrt_deriv3(complex a)
{
    complex r2 = complex_square(complex_inv(a));
    return vec2_scale(complex_mul(r2, complex_inv(complex_sqrt(a))), 0.375);
}

static inline complex
complex_sqrt_deriv4(complex a)
{
    complex r = complex_inv(a);
    complex r3 = complex_mul(complex_square(r), r);
    return vec2_scale(complex_mul(r3, complex_inv(complex_sqrt(a))), -0.9375);
}

static inline complex
complex_sin_deriv1(complex a)
{
    return complex_cos(a);
}

static inline complex
complex_sin_deriv2(complex a)
{
    return complex_neg(complex_sin(a));
}

static inline complex
complex_sin_deriv3(complex a)
{
    return complex_neg(complex_cos(a));
}

static inline complex
complex_sin_deriv4(complex a)
{
    return complex_sin(a);
}

static inline complex
complex_cos_deriv1(complex a)
{
    return complex_neg(complex_sin(a));
}

static inline complex
complex_cos_deriv2(complex a)
{
    return complex_neg(complex_cos(a));
}

static inline complex
complex_cos_deriv3(complex a)
{
    return complex_sin(a);
}

static inline complex
complex_cos_deriv4(complex a)
{
    return complex_cos(a);
}

static inline complex
complex_tan_deriv1(complex a)
{
    return complex_square(complex_inv(complex_cos(a)));
}

static inline complex
complex_tan_deriv2(complex a)
{
    return vec2_scale(complex_mul(complex_tan(a), complex_tan_deriv1(a)), 2);
}

static inline complex
complex_tan_deriv3(complex a)
{
    complex t2 = complex_square(complex_tan(a));
    complex p = complex_radd(1, vec2_scale(t2, 3));
    return vec2_scale(complex_mul(p, complex_tan_deriv1(a)), 2);
}

static inline complex
complex_tan_deriv4(complex a)
{
    complex t = complex_tan(a);
    complex p = complex_radd(2, vec2_scale(complex_square(t), 3));
    return vec2_scale(complex_mul(complex_mul(t, p), complex_tan_deriv1(a)), 8);
}

static inline complex
complex_asin_deriv1(complex a)
{
    /* Separate square roots preserve the principal branches on both cuts. */
    complex p = complex_sqrt(complex_make(1-a.data[0], -a.data[1]));
    complex q = complex_sqrt(complex_make(1+a.data[0], a.data[1]));
    return complex_inv(complex_mul(p, q));
}

static inline complex
complex_asin_deriv2(complex a)
{
    complex q = complex_div(a, complex_rsub(1, complex_square(a)));
    return complex_mul(q, complex_asin_deriv1(a));
}

static inline complex
complex_asin_deriv3(complex a)
{
    complex r = complex_inv(complex_rsub(1, complex_square(a)));
    complex q = complex_mul(a, r);
    complex p = complex_add(complex_square(r), vec2_scale(complex_square(q), 2));
    return complex_mul(p, complex_asin_deriv1(a));
}

static inline complex
complex_asin_deriv4(complex a)
{
    complex r = complex_inv(complex_rsub(1, complex_square(a)));
    complex q = complex_mul(a, r);
    complex p = complex_add(vec2_scale(complex_square(r), 3), vec2_scale(complex_square(q), 2));
    return vec2_scale(complex_mul(complex_mul(q, p), complex_asin_deriv1(a)), 3);
}

static inline complex
complex_acos_deriv1(complex a)
{
    return complex_neg(complex_asin_deriv1(a));
}

static inline complex
complex_acos_deriv2(complex a)
{
    return complex_neg(complex_asin_deriv2(a));
}

static inline complex
complex_acos_deriv3(complex a)
{
    return complex_neg(complex_asin_deriv3(a));
}

static inline complex
complex_acos_deriv4(complex a)
{
    return complex_neg(complex_asin_deriv4(a));
}

static inline complex
complex_atan_deriv1(complex a)
{
    return complex_inv(complex_radd(1, complex_square(a)));
}

static inline complex
complex_atan_deriv2(complex a)
{
    complex r = complex_atan_deriv1(a);
    return vec2_scale(complex_mul(a, complex_square(r)), -2);
}

static inline complex
complex_atan_deriv3(complex a)
{
    complex r = complex_atan_deriv1(a);
    complex q = complex_mul(a, r);
    complex p = complex_sub(vec2_scale(complex_square(q), 3), complex_square(r));
    return vec2_scale(complex_mul(p, r), 2);
}

static inline complex
complex_atan_deriv4(complex a)
{
    complex r = complex_atan_deriv1(a);
    complex q = complex_mul(a, r);
    complex p = complex_sub(complex_square(r), complex_square(q));
    return vec2_scale(complex_mul(complex_mul(q, p), r), 24);
}

static inline complex
complex_sinh_deriv1(complex a)
{
    return complex_cosh(a);
}

static inline complex
complex_sinh_deriv2(complex a)
{
    return complex_sinh(a);
}

static inline complex
complex_sinh_deriv3(complex a)
{
    return complex_cosh(a);
}

static inline complex
complex_sinh_deriv4(complex a)
{
    return complex_sinh(a);
}

static inline complex
complex_cosh_deriv1(complex a)
{
    return complex_sinh(a);
}

static inline complex
complex_cosh_deriv2(complex a)
{
    return complex_cosh(a);
}

static inline complex
complex_cosh_deriv3(complex a)
{
    return complex_sinh(a);
}

static inline complex
complex_cosh_deriv4(complex a)
{
    return complex_cosh(a);
}

static inline complex
complex_tanh_deriv1(complex a)
{
    return complex_square(complex_inv(complex_cosh(a)));
}

static inline complex
complex_tanh_deriv2(complex a)
{
    return vec2_scale(complex_mul(complex_tanh(a), complex_tanh_deriv1(a)), -2);
}

static inline complex
complex_tanh_deriv3(complex a)
{
    complex t2 = complex_square(complex_tanh(a));
    complex p = complex_lsub(vec2_scale(t2, 3), 1);
    return vec2_scale(complex_mul(p, complex_tanh_deriv1(a)), 2);
}

static inline complex
complex_tanh_deriv4(complex a)
{
    complex t = complex_tanh(a);
    complex p = complex_rsub(2, vec2_scale(complex_square(t), 3));
    return vec2_scale(complex_mul(complex_mul(t, p), complex_tanh_deriv1(a)), 8);
}

static inline complex
complex_asinh_deriv1(complex a)
{
    /* sqrt(1+i*a)*sqrt(1-i*a), with signed zeros kept on the cuts. */
    complex p = complex_sqrt(complex_make(1-a.data[1], a.data[0]));
    complex q = complex_sqrt(complex_make(1+a.data[1], -a.data[0]));
    return complex_inv(complex_mul(p, q));
}

static inline complex
complex_asinh_deriv2(complex a)
{
    complex q = complex_div(a, complex_radd(1, complex_square(a)));
    return complex_neg(complex_mul(q, complex_asinh_deriv1(a)));
}

static inline complex
complex_asinh_deriv3(complex a)
{
    complex r = complex_inv(complex_radd(1, complex_square(a)));
    complex q = complex_mul(a, r);
    complex p = complex_sub(vec2_scale(complex_square(q), 2), complex_square(r));
    return complex_mul(p, complex_asinh_deriv1(a));
}

static inline complex
complex_asinh_deriv4(complex a)
{
    complex r = complex_inv(complex_radd(1, complex_square(a)));
    complex q = complex_mul(a, r);
    complex p = complex_sub(vec2_scale(complex_square(r), 3), vec2_scale(complex_square(q), 2));
    return vec2_scale(complex_mul(complex_mul(q, p), complex_asinh_deriv1(a)), 3);
}

static inline complex
complex_acosh_deriv1(complex a)
{
    complex p = complex_sqrt(complex_make(a.data[0]-1, a.data[1]));
    complex q = complex_sqrt(complex_make(a.data[0]+1, a.data[1]));
    return complex_inv(complex_mul(p, q));
}

static inline complex
complex_acosh_deriv2(complex a)
{
    complex q = complex_div(a, complex_lsub(complex_square(a), 1));
    return complex_neg(complex_mul(q, complex_acosh_deriv1(a)));
}

static inline complex
complex_acosh_deriv3(complex a)
{
    complex r = complex_inv(complex_lsub(complex_square(a), 1));
    complex q = complex_mul(a, r);
    complex p = complex_add(vec2_scale(complex_square(q), 2), complex_square(r));
    return complex_mul(p, complex_acosh_deriv1(a));
}

static inline complex
complex_acosh_deriv4(complex a)
{
    complex r = complex_inv(complex_lsub(complex_square(a), 1));
    complex q = complex_mul(a, r);
    complex p = complex_add(vec2_scale(complex_square(q), 2), vec2_scale(complex_square(r), 3));
    return vec2_scale(complex_mul(complex_mul(q, p), complex_acosh_deriv1(a)), -3);
}

static inline complex
complex_atanh_deriv1(complex a)
{
    return complex_inv(complex_rsub(1, complex_square(a)));
}

static inline complex
complex_atanh_deriv2(complex a)
{
    complex r = complex_atanh_deriv1(a);
    return vec2_scale(complex_mul(a, complex_square(r)), 2);
}

static inline complex
complex_atanh_deriv3(complex a)
{
    complex r = complex_atanh_deriv1(a);
    complex q = complex_mul(a, r);
    complex p = complex_add(complex_square(r), vec2_scale(complex_square(q), 3));
    return vec2_scale(complex_mul(p, r), 2);
}

static inline complex
complex_atanh_deriv4(complex a)
{
    complex r = complex_atanh_deriv1(a);
    complex q = complex_mul(a, r);
    complex p = complex_add(complex_square(r), complex_square(q));
    return vec2_scale(complex_mul(complex_mul(q, p), r), 24);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
complex_lmat(complex a, double result[ALGPY_COMPLEX_N*ALGPY_COMPLEX_N])
{
    complex coeffs[ALGPY_COMPLEX_N * ALGPY_COMPLEX_N];
    int row, col, i;
    complex_coeffs(coeffs);
    for (row = 0; row < ALGPY_COMPLEX_N; row++) {
        for (col = 0; col < ALGPY_COMPLEX_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_COMPLEX_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_COMPLEX_N + col].data[row];
            }
            result[row * ALGPY_COMPLEX_N + col] = value;
        }
    }
}

static inline void
complex_rmat(complex a, double result[ALGPY_COMPLEX_N*ALGPY_COMPLEX_N])
{
    complex coeffs[ALGPY_COMPLEX_N * ALGPY_COMPLEX_N];
    int row, col, i;
    complex_coeffs(coeffs);
    for (row = 0; row < ALGPY_COMPLEX_N; row++) {
        for (col = 0; col < ALGPY_COMPLEX_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_COMPLEX_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_COMPLEX_N + i].data[col];
            }
            result[row * ALGPY_COMPLEX_N + col] = value;
        }
    }
}


#endif /* ALGPY_COMPLEX_H */
