/*
 * Core arithmetic of the multicomplex numbers
 * w + x*e1 + y*e2 + z*e1*e2 with e1*e1 = k[0], e2*e2 = k[1].
 *
 * Pure C with no dependency on Python or NumPy. The parameters k are passed
 * explicitly to every operation that depends on them.
 */
#ifndef ALGPY_MULTICOMPLEX2_H
#define ALGPY_MULTICOMPLEX2_H

#include "../../../fields/complex.h"

#define ALGPY_MULTICOMPLEX2_N 4

#include <math.h>

#include "../../multicomplex/complex2.h"
#include "../../multicomplex/splitcomplex2.h"
#include "../../multicomplex/dual2.h"

#include "multicomplex1.h"
#include "../../../vectorspaces/vec1.h"
#include "../../../vectorspaces/vec2.h"
#include "../../../vectorspaces/vec4.h"

typedef vec4 multicomplex2;

static inline multicomplex2
multicomplex2_make(double x, double x0, double x1, double x01)
{
    multicomplex2 r = {{x, x0, x1, x01}};
    return r;
}

static inline multicomplex2
multicomplex2_from_re(double a)
{
    return multicomplex2_make(a, 0, 0, 0);
}

static inline double
multicomplex2_to_re(multicomplex2 a)
{
    return a.data[0];
}

static inline multicomplex2
multicomplex2_one(void)
{
    return multicomplex2_from_re(1);
}

static inline multicomplex2
multicomplex2_mul(vec2 k, multicomplex2 x, multicomplex2 y)
{
    const double *a = x.data, *b = y.data;
    if (vec4_is_real(y)) {
        return vec4_scale(x, b[0]);
    }
    if (vec4_is_real(x)) {
        return vec4_scale(y, a[0]);
    }
    return multicomplex2_make(
        a[0]*b[0] + k.data[0]*a[1]*b[1] + k.data[1]*a[2]*b[2] + k.data[0]*k.data[1]*a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + k.data[1]*a[2]*b[3] + k.data[1]*a[3]*b[2],
        a[0]*b[2] + k.data[0]*a[1]*b[3] + a[2]*b[0] + k.data[0]*a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0]);
}

/* Conjugation of the outer e2 part: (u + v*e2) -> (u - v*e2) */
static inline multicomplex2
multicomplex2_conjugate(multicomplex2 a)
{
    return multicomplex2_make(a.data[0], a.data[1], -a.data[2], -a.data[3]);
}

/* Norm a*conj(a) = u^2 - k[1]*v^2 with a = u + v*e2. */
static inline multicomplex1
multicomplex2_norm2(vec2 k, multicomplex2 a)
{
    vec1 k0 = {{k.data[0]}};
    multicomplex1 u = multicomplex1_make(a.data[0], a.data[1]);
    multicomplex1 v = multicomplex1_make(a.data[2], a.data[3]);
    multicomplex1 uu = multicomplex1_mul(k0, u, u);
    multicomplex1 vv = multicomplex1_mul(k0, v, v);
    return vec2_sub(uu, vec2_scale(vv, k.data[1]));
}

static inline double
multicomplex2_norm4(vec2 k, multicomplex2 a)
{
    vec1 k0 = {{k.data[0]}};
    return multicomplex1_norm2(k0, multicomplex2_norm2(k, a));
}

static inline complex
multicomplex2_norm(vec2 k, multicomplex2 a)
{
    return complex_root_real(multicomplex2_norm4(k, a), 4);
}

static inline multicomplex2
multicomplex2_inv(vec2 k, multicomplex2 a)
{
    vec1 k0 = {{k.data[0]}};
    multicomplex1 ninv = multicomplex1_inv(k0, multicomplex2_norm2(k, a));
    multicomplex1 u = multicomplex1_make(a.data[0], a.data[1]);
    multicomplex1 v = multicomplex1_make(a.data[2], a.data[3]);
    multicomplex1 ru = multicomplex1_mul(k0, u, ninv);
    multicomplex1 rv = multicomplex1_mul(k0, vec2_neg(v), ninv);
    return multicomplex2_make(ru.data[0], ru.data[1], rv.data[0], rv.data[1]);
}

static inline multicomplex2
multicomplex2_div(vec2 k, multicomplex2 a, multicomplex2 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return multicomplex2_mul(k, a, multicomplex2_inv(k, b));
}

static inline void
multicomplex2_basis(multicomplex2 result[ALGPY_MULTICOMPLEX2_N])
{
    result[0] = multicomplex2_make(1, 0, 0, 0);
    result[1] = multicomplex2_make(0, 1, 0, 0);
    result[2] = multicomplex2_make(0, 0, 1, 0);
    result[3] = multicomplex2_make(0, 0, 0, 1);
}

static inline void
multicomplex2_coeffs(
    vec2 k,
    multicomplex2 result[ALGPY_MULTICOMPLEX2_N * ALGPY_MULTICOMPLEX2_N])
{
    multicomplex2 basis[ALGPY_MULTICOMPLEX2_N];
    int i, j;
    multicomplex2_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX2_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX2_N; j++) {
            result[i * ALGPY_MULTICOMPLEX2_N + j] =
                multicomplex2_mul(k, basis[i], basis[j]);
        }
    }
}

/* Exchange the first and outer generators, which commute. */
static inline multicomplex2
multicomplex2_swap_outer(multicomplex2 a)
{
    return multicomplex2_make(a.data[0], a.data[2], a.data[1], a.data[3]);
}

static inline multicomplex2
multicomplex2_apply(vec2 k, multicomplex2 a, multicomplex1 (*f)(vec1, multicomplex1), multicomplex1 (*deriv1)(vec1, multicomplex1))
{
    vec1 base = {{k.data[0]}};
    multicomplex1 u = multicomplex1_make(a.data[0], a.data[1]);
    multicomplex1 v = multicomplex1_make(a.data[2], a.data[3]);
    multicomplex1 p, q;
    multicomplex2 result;
    int i;
    if (k.data[1] == 0) {
        p = f(base, u);
        q = multicomplex1_mul(base, v, deriv1(base, u));
        for (i = 0; i < 2; i++) {
            result.data[i] = p.data[i];
            result.data[i+2] = q.data[i];
        }
        return result;
    }
    if (k.data[1] > 0) {
        double s = sqrt(k.data[1]);
        p = f(base, multicomplex1_add(u, vec2_scale(v, s)));
        q = f(base, multicomplex1_sub(u, vec2_scale(v, s)));
        for (i = 0; i < 2; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+2] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[1] < 0)) {
        return multicomplex2_make(NAN, NAN, NAN, NAN);
    }
    /* Find a complex direction in the lower algebra. */
    int bit = 0;
    double inner = 0;
    if (k.data[0] < 0) {
        bit = 1;
        inner = k.data[0];
    }
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec2 swapped = {{k.data[1], k.data[0]}};
        return multicomplex2_swap_outer(multicomplex2_apply(swapped, multicomplex2_swap_outer(a), f, deriv1));
    }
    double outer_root = sqrt(-k.data[1]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex1 jv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 2; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
    }
    p = f(base, multicomplex1_add(u, jv));
    q = f(base, multicomplex1_sub(u, jv));
    for (i = 0; i < 2; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+2] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex2
multicomplex2_exp(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_exp(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_exp(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_exp(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_exp, multicomplex1_exp_deriv1);
}

static inline multicomplex2
multicomplex2_log(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_log(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_log(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_log(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_log, multicomplex1_log_deriv1);
}

static inline multicomplex2
multicomplex2_pow(vec2 k, multicomplex2 a, multicomplex2 b)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_pow(a, b);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_pow(a, b);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_pow(a, b);
    }
    vec1 base = {{k.data[0]}};
    multicomplex1 u = multicomplex1_make(a.data[0], a.data[1]);
    multicomplex1 v = multicomplex1_make(a.data[2], a.data[3]);
    multicomplex1 bu = multicomplex1_make(b.data[0], b.data[1]);
    multicomplex1 bv = multicomplex1_make(b.data[2], b.data[3]);
    multicomplex1 p, q;
    multicomplex2 result;
    int i;
    if (k.data[1] == 0) {
        return multicomplex2_exp(k, multicomplex2_mul(k, multicomplex2_log(k, a), b));
    }
    if (k.data[1] > 0) {
        double s = sqrt(k.data[1]);
        p = multicomplex1_pow(base, multicomplex1_add(u, vec2_scale(v, s)), multicomplex1_add(bu, vec2_scale(bv, s)));
        q = multicomplex1_pow(base, multicomplex1_sub(u, vec2_scale(v, s)), multicomplex1_sub(bu, vec2_scale(bv, s)));
        for (i = 0; i < 2; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+2] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[1] < 0)) {
        return multicomplex2_make(NAN, NAN, NAN, NAN);
    }
    /* Find a complex direction in the lower algebra. */
    int bit = 0;
    double inner = 0;
    if (k.data[0] < 0) {
        bit = 1;
        inner = k.data[0];
    }
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec2 swapped = {{k.data[1], k.data[0]}};
        return multicomplex2_swap_outer(multicomplex2_pow(swapped, multicomplex2_swap_outer(a), multicomplex2_swap_outer(b)));
    }
    double outer_root = sqrt(-k.data[1]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex1 jv;
    multicomplex1 jbv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 2; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
        jbv.data[i] = bv.data[i^bit]*(i & bit ? odd : even);
    }
    p = multicomplex1_pow(base, multicomplex1_add(u, jv), multicomplex1_add(bu, jbv));
    q = multicomplex1_pow(base, multicomplex1_sub(u, jv), multicomplex1_sub(bu, jbv));
    for (i = 0; i < 2; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+2] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex2
multicomplex2_sqrt(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_sqrt(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_sqrt(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_sqrt(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_sqrt, multicomplex1_sqrt_deriv1);
}

static inline void
multicomplex2_minpoly(vec2 k, multicomplex2 a, double result[5])
{
    double k0 = k.data[0], k1 = k.data[1];
    double x0 = a.data[0], x1 = a.data[1], x2 = a.data[2], x3 = a.data[3];
    double x02 = x0*x0, x12 = x1*x1, x22 = x2*x2, x32 = x3*x3;
    result[4] = 1;
    result[3] = -4*x0;
    result[2] = -2*k0*k1*x32 - 2*k0*x12 - 2*k1*x22 + 6*x02;
    result[1] = 4*k0*k1*x0*x32 - 8*k0*k1*x1*x2*x3 + 4*k0*x0*x12 + 4*k1*x0*x22 - 4*x02*x0;
    result[0] = k0*k0*k1*k1*x32*x32 - 2*k0*k0*k1*x12*x32 + k0*k0*x12*x12
        - 2*k0*k1*k1*x22*x32 - 2*k0*k1*x02*x32 + 8*k0*k1*x0*x1*x2*x3
        - 2*k0*k1*x12*x22 - 2*k0*x02*x12 + k1*k1*x22*x22 - 2*k1*x02*x22
        + x02*x02;
}

/*
 * Eigenvalues x0 +/- sqrt(k0)*x1 +/- sqrt(k1)*x2 +/- sqrt(k0)*sqrt(k1)*x3,
 * as complex numbers.
 */
static inline void
multicomplex2_eivals(vec2 k, multicomplex2 a, complex result[4])
{
    complex s0 = complex_sqrt_real(k.data[0]);
    complex s1 = complex_sqrt_real(k.data[1]);
    complex t1 = complex_mul(complex_make(a.data[1], 0), s0);
    complex t2 = complex_mul(complex_make(a.data[2], 0), s1);
    complex t3 = complex_mul(complex_mul(complex_make(a.data[3], 0), s0), s1);
    result[0] = complex_make(a.data[0] - t1.data[0] - t2.data[0] + t3.data[0], -t1.data[1] - t2.data[1] + t3.data[1]);
    result[1] = complex_make(a.data[0] + t1.data[0] + t2.data[0] + t3.data[0], +t1.data[1] + t2.data[1] + t3.data[1]);
    result[2] = complex_make(a.data[0] + t1.data[0] - t2.data[0] - t3.data[0], +t1.data[1] - t2.data[1] - t3.data[1]);
    result[3] = complex_make(a.data[0] - t1.data[0] + t2.data[0] - t3.data[0], -t1.data[1] + t2.data[1] - t3.data[1]);
}

/* Common algebra API. */
static inline multicomplex2
multicomplex2_zero(void)
{
    return vec4_zero();
}

static inline int
multicomplex2_is_real(multicomplex2 a)
{
    return vec4_is_real(a);
}

static inline multicomplex2
multicomplex2_pos(multicomplex2 a)
{
    return vec4_pos(a);
}

static inline multicomplex2
multicomplex2_neg(multicomplex2 a)
{
    return vec4_neg(a);
}

static inline multicomplex2
multicomplex2_add(multicomplex2 a, multicomplex2 b)
{
    return vec4_add(a, b);
}

static inline multicomplex2
multicomplex2_sub(multicomplex2 a, multicomplex2 b)
{
    return vec4_sub(a, b);
}

static inline int
multicomplex2_equal(multicomplex2 a, multicomplex2 b)
{
    return vec4_all_eq(a, b);
}

static inline int
multicomplex2_not_equal(multicomplex2 a, multicomplex2 b)
{
    return !multicomplex2_equal(a, b);
}

static inline multicomplex2
multicomplex2_re(multicomplex2 a)
{
    return multicomplex2_from_re(a.data[0]);
}

static inline multicomplex2
multicomplex2_im(multicomplex2 a)
{
    return multicomplex2_sub(a, multicomplex2_re(a));
}

static inline multicomplex2
multicomplex2_cos(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_cos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_cos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_cos(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_cos, multicomplex1_cos_deriv1);
}

static inline multicomplex2
multicomplex2_sin(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_sin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_sin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_sin(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_sin, multicomplex1_sin_deriv1);
}

static inline multicomplex2
multicomplex2_tan(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_tan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_tan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_tan(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_tan, multicomplex1_tan_deriv1);
}

static inline multicomplex2
multicomplex2_acos(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_acos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_acos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_acos(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_acos, multicomplex1_acos_deriv1);
}

static inline multicomplex2
multicomplex2_asin(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_asin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_asin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_asin(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_asin, multicomplex1_asin_deriv1);
}

static inline multicomplex2
multicomplex2_atan(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_atan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_atan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_atan(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_atan, multicomplex1_atan_deriv1);
}

static inline multicomplex2
multicomplex2_cosh(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_cosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_cosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_cosh(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_cosh, multicomplex1_cosh_deriv1);
}

static inline multicomplex2
multicomplex2_sinh(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_sinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_sinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_sinh(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_sinh, multicomplex1_sinh_deriv1);
}

static inline multicomplex2
multicomplex2_tanh(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_tanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_tanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_tanh(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_tanh, multicomplex1_tanh_deriv1);
}

static inline multicomplex2
multicomplex2_acosh(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_acosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_acosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_acosh(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_acosh, multicomplex1_acosh_deriv1);
}

static inline multicomplex2
multicomplex2_asinh(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_asinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_asinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_asinh(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_asinh, multicomplex1_asinh_deriv1);
}

static inline multicomplex2
multicomplex2_atanh(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_atanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_atanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_atanh(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_atanh, multicomplex1_atanh_deriv1);
}

static inline multicomplex2
multicomplex2_square(vec2 k, multicomplex2 a)
{
    return multicomplex2_mul(k, a, a);
}

static inline multicomplex2
multicomplex2_radd(double a, multicomplex2 b)
{
    return multicomplex2_add(multicomplex2_from_re(a), b);
}

static inline multicomplex2
multicomplex2_ladd(multicomplex2 a, double b)
{
    return multicomplex2_add(a, multicomplex2_from_re(b));
}

static inline multicomplex2
multicomplex2_rsub(double a, multicomplex2 b)
{
    return multicomplex2_sub(multicomplex2_from_re(a), b);
}

static inline multicomplex2
multicomplex2_lsub(multicomplex2 a, double b)
{
    return multicomplex2_sub(a, multicomplex2_from_re(b));
}

static inline multicomplex2
multicomplex2_rmul(vec2 k, double a, multicomplex2 b)
{
    return multicomplex2_mul(k, multicomplex2_from_re(a), b);
}

static inline multicomplex2
multicomplex2_lmul(vec2 k, multicomplex2 a, double b)
{
    return multicomplex2_mul(k, a, multicomplex2_from_re(b));
}

static inline multicomplex2
multicomplex2_rdiv(vec2 k, double a, multicomplex2 b)
{
    return multicomplex2_div(k, multicomplex2_from_re(a), b);
}

static inline multicomplex2
multicomplex2_ldiv(vec2 k, multicomplex2 a, double b)
{
    return multicomplex2_div(k, a, multicomplex2_from_re(b));
}

static inline multicomplex2
multicomplex2_rpow(vec2 k, double a, multicomplex2 b)
{
    return multicomplex2_pow(k, multicomplex2_from_re(a), b);
}

static inline multicomplex2
multicomplex2_lpow(vec2 k, multicomplex2 a, double b)
{
    return multicomplex2_pow(k, a, multicomplex2_from_re(b));
}

static inline multicomplex2
multicomplex2_expm1(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_expm1(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_expm1(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_expm1(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_expm1, multicomplex1_expm1_deriv1);
}

static inline multicomplex2
multicomplex2_log1p(vec2 k, multicomplex2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return complex2_log1p(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1) {
        return splitcomplex2_log1p(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0) {
        return dual2_log1p(a);
    }
    return multicomplex2_apply(k, a, multicomplex1_log1p, multicomplex1_log1p_deriv1);
}

static inline multicomplex2
multicomplex2_mix(vec2 k, multicomplex2 a, multicomplex2 b, double t)
{
    multicomplex2 delta = multicomplex2_mul(k, multicomplex2_inv(k, a), b);
    return multicomplex2_mul(k, a, multicomplex2_lpow(k, delta, t));
}

/* Derivatives through order 2, for nilpotent extensions at higher levels. */
static inline multicomplex2
multicomplex2_exp_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_exp_deriv1, multicomplex1_exp_deriv2);
}

static inline multicomplex2
multicomplex2_exp_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_exp_deriv2, multicomplex1_exp_deriv3);
}

static inline multicomplex2
multicomplex2_expm1_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_expm1_deriv1, multicomplex1_expm1_deriv2);
}

static inline multicomplex2
multicomplex2_expm1_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_expm1_deriv2, multicomplex1_expm1_deriv3);
}

static inline multicomplex2
multicomplex2_log_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_log_deriv1, multicomplex1_log_deriv2);
}

static inline multicomplex2
multicomplex2_log_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_log_deriv2, multicomplex1_log_deriv3);
}

static inline multicomplex2
multicomplex2_log1p_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_log1p_deriv1, multicomplex1_log1p_deriv2);
}

static inline multicomplex2
multicomplex2_log1p_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_log1p_deriv2, multicomplex1_log1p_deriv3);
}

static inline multicomplex2
multicomplex2_inv_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_inv_deriv1, multicomplex1_inv_deriv2);
}

static inline multicomplex2
multicomplex2_inv_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_inv_deriv2, multicomplex1_inv_deriv3);
}

static inline multicomplex2
multicomplex2_square_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_square_deriv1, multicomplex1_square_deriv2);
}

static inline multicomplex2
multicomplex2_square_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_square_deriv2, multicomplex1_square_deriv3);
}

static inline multicomplex2
multicomplex2_sqrt_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_sqrt_deriv1, multicomplex1_sqrt_deriv2);
}

static inline multicomplex2
multicomplex2_sqrt_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_sqrt_deriv2, multicomplex1_sqrt_deriv3);
}

static inline multicomplex2
multicomplex2_sin_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_sin_deriv1, multicomplex1_sin_deriv2);
}

static inline multicomplex2
multicomplex2_sin_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_sin_deriv2, multicomplex1_sin_deriv3);
}

static inline multicomplex2
multicomplex2_cos_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_cos_deriv1, multicomplex1_cos_deriv2);
}

static inline multicomplex2
multicomplex2_cos_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_cos_deriv2, multicomplex1_cos_deriv3);
}

static inline multicomplex2
multicomplex2_tan_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_tan_deriv1, multicomplex1_tan_deriv2);
}

static inline multicomplex2
multicomplex2_tan_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_tan_deriv2, multicomplex1_tan_deriv3);
}

static inline multicomplex2
multicomplex2_asin_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_asin_deriv1, multicomplex1_asin_deriv2);
}

static inline multicomplex2
multicomplex2_asin_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_asin_deriv2, multicomplex1_asin_deriv3);
}

static inline multicomplex2
multicomplex2_acos_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_acos_deriv1, multicomplex1_acos_deriv2);
}

static inline multicomplex2
multicomplex2_acos_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_acos_deriv2, multicomplex1_acos_deriv3);
}

static inline multicomplex2
multicomplex2_atan_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_atan_deriv1, multicomplex1_atan_deriv2);
}

static inline multicomplex2
multicomplex2_atan_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_atan_deriv2, multicomplex1_atan_deriv3);
}

static inline multicomplex2
multicomplex2_sinh_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_sinh_deriv1, multicomplex1_sinh_deriv2);
}

static inline multicomplex2
multicomplex2_sinh_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_sinh_deriv2, multicomplex1_sinh_deriv3);
}

static inline multicomplex2
multicomplex2_cosh_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_cosh_deriv1, multicomplex1_cosh_deriv2);
}

static inline multicomplex2
multicomplex2_cosh_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_cosh_deriv2, multicomplex1_cosh_deriv3);
}

static inline multicomplex2
multicomplex2_tanh_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_tanh_deriv1, multicomplex1_tanh_deriv2);
}

static inline multicomplex2
multicomplex2_tanh_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_tanh_deriv2, multicomplex1_tanh_deriv3);
}

static inline multicomplex2
multicomplex2_asinh_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_asinh_deriv1, multicomplex1_asinh_deriv2);
}

static inline multicomplex2
multicomplex2_asinh_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_asinh_deriv2, multicomplex1_asinh_deriv3);
}

static inline multicomplex2
multicomplex2_acosh_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_acosh_deriv1, multicomplex1_acosh_deriv2);
}

static inline multicomplex2
multicomplex2_acosh_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_acosh_deriv2, multicomplex1_acosh_deriv3);
}

static inline multicomplex2
multicomplex2_atanh_deriv1(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_atanh_deriv1, multicomplex1_atanh_deriv2);
}

static inline multicomplex2
multicomplex2_atanh_deriv2(vec2 k, multicomplex2 a)
{
    return multicomplex2_apply(k, a, multicomplex1_atanh_deriv2, multicomplex1_atanh_deriv3);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
multicomplex2_lmat(vec2 k, multicomplex2 a, double result[ALGPY_MULTICOMPLEX2_N*ALGPY_MULTICOMPLEX2_N])
{
    multicomplex2 coeffs[ALGPY_MULTICOMPLEX2_N * ALGPY_MULTICOMPLEX2_N];
    int row, col, i;
    multicomplex2_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_MULTICOMPLEX2_N; row++) {
        for (col = 0; col < ALGPY_MULTICOMPLEX2_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_MULTICOMPLEX2_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_MULTICOMPLEX2_N + col].data[row];
            }
            result[row * ALGPY_MULTICOMPLEX2_N + col] = value;
        }
    }
}

static inline void
multicomplex2_rmat(vec2 k, multicomplex2 a, double result[ALGPY_MULTICOMPLEX2_N*ALGPY_MULTICOMPLEX2_N])
{
    multicomplex2 coeffs[ALGPY_MULTICOMPLEX2_N * ALGPY_MULTICOMPLEX2_N];
    int row, col, i;
    multicomplex2_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_MULTICOMPLEX2_N; row++) {
        for (col = 0; col < ALGPY_MULTICOMPLEX2_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_MULTICOMPLEX2_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_MULTICOMPLEX2_N + i].data[col];
            }
            result[row * ALGPY_MULTICOMPLEX2_N + col] = value;
        }
    }
}


#endif /* ALGPY_MULTICOMPLEX2_H */
