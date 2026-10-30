/*
 * Shared arithmetic of complex numbers (x + y*i), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX_H
#define ALGPY_COMPLEX_H

#include <stdint.h>
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

/* MSVC's complex math API uses a struct rather than the C99 complex type. */
#ifdef _MSC_VER
typedef _Dcomplex algpy_native_complex;
#else
typedef double _Complex algpy_native_complex;
#endif

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
    double x = a.data[0], y = a.data[1];
    /* Normalize by the larger component, keeping the squared ratio <= 1.
     * Divide by that component last so intermediate reciprocals stay finite. */
    if (fabs(x) >= fabs(y)) {
        double ratio = y/x;
        double scale = (1/(1 + ratio*ratio))/x;
        return complex_make(scale, -ratio*scale);
    }
    double ratio = x/y;
    double scale = (1/(1 + ratio*ratio))/y;
    return complex_make(ratio*scale, -scale);
}

static inline complex
complex_div(complex a, complex b)
{
    double ar = a.data[0], ai = a.data[1], br = b.data[0], bi = b.data[1];
    double amax = fmax(fabs(ar), fabs(ai)), bmax = fmax(fabs(br), fabs(bi));
    double denominator;
    int ae = 0, be = 0;
    if (bi == 0 && br != 0) {
        return complex_make(ar/br, ai/br);
    }
    if (br == 0 && bi != 0) {
        return complex_make(ai/bi, -ar/bi);
    }
    if (!isfinite(amax) || !isfinite(bmax) || bmax == 0) {
        return complex_mul(a, complex_inv(b));
    }
    /* Scale both operands by exact powers of two before forming products.
     * This also covers subnormal divisors whose inverse is unrepresentable. */
    frexp(amax, &ae);
    frexp(bmax, &be);
    ar = scalbn(ar, -ae); ai = scalbn(ai, -ae);
    br = scalbn(br, -be); bi = scalbn(bi, -be);
    denominator = br*br + bi*bi;
    return complex_make(scalbn((ar*br + ai*bi)/denominator, ae - be),
                        scalbn((ai*br - ar*bi)/denominator, ae - be));
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

/* Integer powers by squaring; negative exponents invert the base once. */
static inline complex
complex_ipow(complex a, int64_t b)
{
    complex result = complex_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return complex_mul(a, a);
    }
    if (b == -1) {
        return complex_inv(a);
    }
    if (b < 0) {
        a = complex_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = complex_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = complex_mul(a, a);
        }
    }
    return result;
}

static inline complex
complex_rpow(complex a, complex b)
{
    return complex_exp(complex_mul(complex_log(a), b));
}

static inline complex
complex_pow(complex a, complex b)
{
    return complex_rpow(a, b);
}

static inline complex
complex_lpow(complex a, complex b)
{
    return complex_rpow(a, b);
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
complex_rpows(double a, complex b) {
  return complex_rpow(complex_from_re(a), b);
}

static inline complex
complex_lpows(complex a, double b) {
  return complex_rpow(a, complex_from_re(b));
}

static inline complex
complex_mix(complex a, complex b, double t)
{
    complex delta = complex_mul(complex_inv(a), b);
    return complex_mul(a, complex_lpows(delta, t));
}

static inline complex
complex_expm1(complex a)
{
    double half_sin = sin(0.5*a.data[1]);
    return complex_make(expm1(a.data[0])*cos(a.data[1]) - 2*half_sin*half_sin,
                        exp(a.data[0])*sin(a.data[1]));
}

static inline complex
complex_log1p(complex a)
{
    if (fabs(a.data[0]) < 0.5 && fabs(a.data[1]) < 0.5) {
        double real = a.data[1] == 0 ? log1p(a.data[0]) :
            0.5*log1p(a.data[0]*(2 + a.data[0]) + a.data[1]*a.data[1]);
        return complex_make(real, atan2(a.data[1], 1 + a.data[0]));
    }
    return complex_log(complex_ladd(a, 1));
}

static inline complex
complex_sqrt(complex a)
{
    algpy_native_complex z;
    /* The native implementation avoids cancellation near the real axis. */
    memcpy(&z, a.data, sizeof(z));
    z = csqrt(z);
    return complex_make(creal(z), cimag(z));
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
    algpy_native_complex z;
    /* Native reduction avoids overflowing sin/cos for large imaginary parts. */
    memcpy(&z, a.data, sizeof(z));
    z = ctan(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_acos(complex a)
{
    algpy_native_complex z;
    /* Both native types store consecutive real and imaginary doubles.
     * Copying preserves signed zeros on branch cuts without arithmetic. */
    memcpy(&z, a.data, sizeof(z));
    z = cacos(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_asin(complex a)
{
    algpy_native_complex z;
    memcpy(&z, a.data, sizeof(z));
    z = casin(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_atan(complex a)
{
    algpy_native_complex z;
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
    algpy_native_complex z;
    /* The bounded result remains finite when sinh/cosh would overflow. */
    memcpy(&z, a.data, sizeof(z));
    z = ctanh(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_acosh(complex a)
{
    algpy_native_complex z;
    memcpy(&z, a.data, sizeof(z));
    z = cacosh(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_asinh(complex a)
{
    algpy_native_complex z;
    memcpy(&z, a.data, sizeof(z));
    z = casinh(z);
    return complex_make(creal(z), cimag(z));
}

static inline complex
complex_atanh(complex a)
{
    algpy_native_complex z;
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
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
complex_lmat(complex a, double result[ALGPY_COMPLEX_N*ALGPY_COMPLEX_N])
{
    complex basis = {{0}};
    complex product;
    int row, col;
    for (col = 0; col < ALGPY_COMPLEX_N; col++) {
        basis.data[col] = 1;
        product = complex_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_COMPLEX_N; row++) {
            result[row * ALGPY_COMPLEX_N + col] = product.data[row];
        }
    }
}

static inline void
complex_rmat(complex a, double result[ALGPY_COMPLEX_N*ALGPY_COMPLEX_N])
{
    complex basis = {{0}};
    complex product;
    int row, col;
    for (row = 0; row < ALGPY_COMPLEX_N; row++) {
        basis.data[row] = 1;
        product = complex_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_COMPLEX_N; col++) {
            result[row * ALGPY_COMPLEX_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_COMPLEX_H */
