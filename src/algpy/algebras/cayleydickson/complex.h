/*
 * Shared arithmetic of complex numbers (x + y*i), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX_H
#define ALGPY_COMPLEX_H

#include <math.h>

#include "../../vectorspaces/vec2.h"

#define ALGPY_COMPLEX_N ALGPY_VEC2_N

/* C99 <complex.h> (pulled in by NumPy) defines complex as a macro. */
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
    result[0] = complex_make(1, 0);
    result[1] = complex_make(0, 1);
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
complex_sqrt(complex a)
{
    double r = hypot(a.data[0], a.data[1]);
    double re = sqrt(fmax(0, (r + a.data[0]) / 2));
    double im = sqrt(fmax(0, (r - a.data[0]) / 2));
    return complex_make(re, copysign(im, a.data[1]));
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

#endif /* ALGPY_COMPLEX_H */