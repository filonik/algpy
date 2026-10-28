/*
 * Core arithmetic of the multicomplex numbers
 * sum_S x_S * e_S over subsets S of {0, 1, 2, 3} (component index = bitmask of S)
 * with e0*e0 = k[0], e1*e1 = k[1], e2*e2 = k[2], e3*e3 = k[3].
 *
 * Pure C with no dependency on Python or NumPy. The parameters k are passed
 * explicitly to every operation that depends on them.
 */
#ifndef ALGPY_MULTICOMPLEX4_H
#define ALGPY_MULTICOMPLEX4_H

#include "../../../fields/complex.h"

#define ALGPY_MULTICOMPLEX4_N 16

#include <math.h>

#include "../../multicomplex/complex4.h"
#include "../../multicomplex/splitcomplex4.h"
#include "../../multicomplex/dual4.h"

#include "multicomplex3.h"
#include "../../../vectorspaces/vec4.h"
#include "../../../vectorspaces/vec16.h"

typedef vec16 multicomplex4;

static inline multicomplex4
multicomplex4_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012,
                    double x3, double x03, double x13, double x013,
                    double x23, double x023, double x123, double x0123)
{
    multicomplex4 r = {{x, x0, x1, x01, x2, x02, x12, x012,
                         x3, x03, x13, x013, x23, x023, x123, x0123}};
    return r;
}

static inline multicomplex4
multicomplex4_from_re(double a)
{
    return multicomplex4_make(a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
multicomplex4_to_re(multicomplex4 a)
{
    return a.data[0];
}

static inline multicomplex4
multicomplex4_one(void)
{
    return multicomplex4_from_re(1);
}

static inline multicomplex4
multicomplex4_rsub(double a, multicomplex4 b) {
  return vec16_sub(multicomplex4_from_re(a), b);
}

static inline multicomplex4
multicomplex4_lsub(multicomplex4 a, double b) {
  return vec16_sub(a, multicomplex4_from_re(b));
}

/* e_i * e_j = (product of k[b] over the bits b shared by i and j) * e_(i^j). */
static inline multicomplex4
multicomplex4_mul(vec4 h, multicomplex4 a, multicomplex4 b)
{
    const double *k = h.data;
    multicomplex4 r;
    int i, j, c;
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        r.data[i] = 0;
    }
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX4_N; j++) {
            double t = a.data[i] * b.data[j];
            for (c = 0; c < 4; c++) {
                if ((i & j) & (1 << c)) {
                    t *= k[c];
                }
            }
            r.data[i ^ j] += t;
        }
    }
    return r;
}

static inline multicomplex4
multicomplex4_conjugate(multicomplex4 a)
{
    multicomplex4 r = a;
    int i;
    for (i = 8; i < ALGPY_MULTICOMPLEX4_N; i++) {
        r.data[i] = -r.data[i];
    }
    return r;
}

static inline multicomplex3
multicomplex4_lower(multicomplex4 a)
{
    return multicomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3],
                              a.data[4], a.data[5], a.data[6], a.data[7]);
}

static inline multicomplex3
multicomplex4_upper(multicomplex4 a)
{
    return multicomplex3_make(a.data[8], a.data[9], a.data[10], a.data[11],
                              a.data[12], a.data[13], a.data[14], a.data[15]);
}

static inline multicomplex4
multicomplex4_join(multicomplex3 u, multicomplex3 v)
{
    return multicomplex4_make(u.data[0], u.data[1], u.data[2], u.data[3],
                              u.data[4], u.data[5], u.data[6], u.data[7],
                              v.data[0], v.data[1], v.data[2], v.data[3],
                              v.data[4], v.data[5], v.data[6], v.data[7]);
}

/* Norm a*conj(a) = u^2 - k[3]*v^2 with a = u + v*e3. */
static inline multicomplex3
multicomplex4_norm2(vec4 k, multicomplex4 a)
{
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    multicomplex3 u = multicomplex4_lower(a);
    multicomplex3 v = multicomplex4_upper(a);
    multicomplex3 uu = multicomplex3_mul(k012, u, u);
    multicomplex3 vv = multicomplex3_mul(k012, v, v);
    return vec8_sub(uu, vec8_scale(vv, k.data[3]));
}

static inline multicomplex2
multicomplex4_norm4(vec4 k, multicomplex4 a)
{
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    return multicomplex3_norm2(k012, multicomplex4_norm2(k, a));
}

static inline multicomplex1
multicomplex4_norm8(vec4 k, multicomplex4 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    return multicomplex2_norm2(k01, multicomplex4_norm4(k, a));
}

static inline double
multicomplex4_norm16(vec4 k, multicomplex4 a)
{
    vec1 k0 = {{k.data[0]}};
    return multicomplex1_norm2(k0, multicomplex4_norm8(k, a));
}

static inline complex
multicomplex4_norm(vec4 k, multicomplex4 a)
{
    return complex_root_real(multicomplex4_norm16(k, a), 16);
}

static inline multicomplex4
multicomplex4_inv(vec4 k, multicomplex4 a)
{
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    multicomplex3 ninv = multicomplex3_inv(k012, multicomplex4_norm2(k, a));
    multicomplex3 u = multicomplex4_lower(a);
    multicomplex3 v = multicomplex4_upper(a);
    multicomplex3 ru = multicomplex3_mul(k012, u, ninv);
    multicomplex3 rv = multicomplex3_mul(k012, vec8_neg(v), ninv);
    return multicomplex4_join(ru, rv);
}

static inline multicomplex4
multicomplex4_div(vec4 k, multicomplex4 a, multicomplex4 b)
{
    if (vec16_is_real(b)) {
        return vec16_div_scalar(a, b.data[0]);
    }
    return multicomplex4_mul(k, a, multicomplex4_inv(k, b));
}

static inline void
multicomplex4_basis(multicomplex4 result[ALGPY_MULTICOMPLEX4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
multicomplex4_coeffs(
    vec4 k,
    multicomplex4 result[ALGPY_MULTICOMPLEX4_N * ALGPY_MULTICOMPLEX4_N])
{
    multicomplex4 basis[ALGPY_MULTICOMPLEX4_N];
    int i, j;
    multicomplex4_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX4_N; j++) {
            result[i * ALGPY_MULTICOMPLEX4_N + j] =
                multicomplex4_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
multicomplex4_minpoly(vec4 k, multicomplex4 a, double result[17])
{
    const int n = ALGPY_MULTICOMPLEX4_N;
    multicomplex4 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * multicomplex4_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            multicomplex4 adj = multicomplex4_lsub(uk, ck);
            uk = multicomplex4_mul(k, a, adj);
        }
    }
}

/*
 * Eigenvalues sum of x_S * prod_{i in S} (+/- sqrt(k_i)) over all subsets S,
 * with the sign of sqrt(k_i) given by the bits of the eigenvalue index,
 * as complex numbers.
 */
static inline void
multicomplex4_eivals(vec4 k, multicomplex4 a, complex result[16])
{
    complex s[4];
    complex n[16];
    int i, j, b;

    for (b = 0; b < 4; b++) {
        s[b] = complex_sqrt_real(k.data[b]);
    }
    for (i = 0; i < 16; i++) {
        complex t = complex_make(a.data[i], 0);
        for (b = 0; b < 4; b++) {
            if (i & (1 << b)) {
                t = complex_mul(t, s[b]);
            }
        }
        n[i] = t;
    }
    for (j = 0; j < 16; j++) {
        double re = 0, im = 0;
        for (i = 0; i < 16; i++) {
            int sign = 1;
            for (b = 0; b < 4; b++) {
                if ((i & (1 << b)) && !(j & (1 << b))) {
                    sign = -sign;
                }
            }
            re += sign * n[i].data[0];
            im += sign * n[i].data[1];
        }
            result[j] = complex_make(re, im);
    }
}

/* Common algebra API. */
static inline multicomplex4
multicomplex4_zero(void)
{
    return vec16_zero();
}

static inline int
multicomplex4_is_real(multicomplex4 a)
{
    return vec16_is_real(a);
}

static inline multicomplex4
multicomplex4_pos(multicomplex4 a)
{
    return vec16_pos(a);
}

static inline multicomplex4
multicomplex4_neg(multicomplex4 a)
{
    return vec16_neg(a);
}

static inline multicomplex4
multicomplex4_add(multicomplex4 a, multicomplex4 b)
{
    return vec16_add(a, b);
}

static inline multicomplex4
multicomplex4_sub(multicomplex4 a, multicomplex4 b)
{
    return vec16_sub(a, b);
}

static inline int
multicomplex4_equal(multicomplex4 a, multicomplex4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
multicomplex4_not_equal(multicomplex4 a, multicomplex4 b)
{
    return !multicomplex4_equal(a, b);
}

static inline multicomplex4
multicomplex4_re(multicomplex4 a)
{
    return multicomplex4_from_re(a.data[0]);
}

static inline multicomplex4
multicomplex4_im(multicomplex4 a)
{
    return multicomplex4_sub(a, multicomplex4_re(a));
}

/* Exchange the first and outer generators, which commute. */
static inline multicomplex4
multicomplex4_swap_outer(multicomplex4 a)
{
    return multicomplex4_make(a.data[0], a.data[8], a.data[2], a.data[10], a.data[4], a.data[12], a.data[6], a.data[14], a.data[1], a.data[9], a.data[3], a.data[11], a.data[5], a.data[13], a.data[7], a.data[15]);
}

static inline multicomplex4
multicomplex4_apply(vec4 k, multicomplex4 a, multicomplex3 (*f)(vec3, multicomplex3), multicomplex3 (*deriv1)(vec3, multicomplex3))
{
    vec3 base = {{k.data[0], k.data[1], k.data[2]}};
    multicomplex3 u = multicomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    multicomplex3 v = multicomplex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    multicomplex3 p, q;
    multicomplex4 result;
    int i;
    if (k.data[3] == 0) {
        p = f(base, u);
        q = multicomplex3_mul(base, v, deriv1(base, u));
        for (i = 0; i < 8; i++) {
            result.data[i] = p.data[i];
            result.data[i+8] = q.data[i];
        }
        return result;
    }
    if (k.data[3] > 0) {
        double s = sqrt(k.data[3]);
        p = f(base, multicomplex3_add(u, vec8_scale(v, s)));
        q = f(base, multicomplex3_sub(u, vec8_scale(v, s)));
        for (i = 0; i < 8; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+8] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[3] < 0)) {
        return multicomplex4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
    }
    /* Find a complex direction in the lower algebra. */
    int bit = 0;
    double inner = 0;
    if (k.data[0] < 0) {
        bit = 1;
        inner = k.data[0];
    }
    else if (k.data[1] < 0) {
        bit = 2;
        inner = k.data[1];
    }
    else if (k.data[2] < 0) {
        bit = 4;
        inner = k.data[2];
    }
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec4 swapped = {{k.data[3], k.data[1], k.data[2], k.data[0]}};
        return multicomplex4_swap_outer(multicomplex4_apply(swapped, multicomplex4_swap_outer(a), f, deriv1));
    }
    double outer_root = sqrt(-k.data[3]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex3 jv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 8; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
    }
    p = f(base, multicomplex3_add(u, jv));
    q = f(base, multicomplex3_sub(u, jv));
    for (i = 0; i < 8; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+8] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex4
multicomplex4_exp(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_exp(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_exp(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_exp(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_exp, multicomplex3_exp_deriv1);
}

static inline multicomplex4
multicomplex4_log(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_log(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_log(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_log(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_log, multicomplex3_log_deriv1);
}

static inline multicomplex4
multicomplex4_sqrt(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_sqrt(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_sqrt(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_sqrt(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_sqrt, multicomplex3_sqrt_deriv1);
}

static inline multicomplex4
multicomplex4_cos(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_cos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_cos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_cos(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_cos, multicomplex3_cos_deriv1);
}

static inline multicomplex4
multicomplex4_sin(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_sin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_sin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_sin(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_sin, multicomplex3_sin_deriv1);
}

static inline multicomplex4
multicomplex4_tan(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_tan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_tan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_tan(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_tan, multicomplex3_tan_deriv1);
}

static inline multicomplex4
multicomplex4_acos(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_acos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_acos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_acos(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_acos, multicomplex3_acos_deriv1);
}

static inline multicomplex4
multicomplex4_asin(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_asin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_asin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_asin(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_asin, multicomplex3_asin_deriv1);
}

static inline multicomplex4
multicomplex4_atan(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_atan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_atan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_atan(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_atan, multicomplex3_atan_deriv1);
}

static inline multicomplex4
multicomplex4_cosh(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_cosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_cosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_cosh(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_cosh, multicomplex3_cosh_deriv1);
}

static inline multicomplex4
multicomplex4_sinh(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_sinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_sinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_sinh(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_sinh, multicomplex3_sinh_deriv1);
}

static inline multicomplex4
multicomplex4_tanh(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_tanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_tanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_tanh(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_tanh, multicomplex3_tanh_deriv1);
}

static inline multicomplex4
multicomplex4_acosh(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_acosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_acosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_acosh(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_acosh, multicomplex3_acosh_deriv1);
}

static inline multicomplex4
multicomplex4_asinh(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_asinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_asinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_asinh(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_asinh, multicomplex3_asinh_deriv1);
}

static inline multicomplex4
multicomplex4_atanh(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_atanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_atanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_atanh(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_atanh, multicomplex3_atanh_deriv1);
}

static inline multicomplex4
multicomplex4_pow(vec4 k, multicomplex4 a, multicomplex4 b)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_pow(a, b);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_pow(a, b);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_pow(a, b);
    }
    vec3 base = {{k.data[0], k.data[1], k.data[2]}};
    multicomplex3 u = multicomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    multicomplex3 v = multicomplex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    multicomplex3 bu = multicomplex3_make(b.data[0], b.data[1], b.data[2], b.data[3], b.data[4], b.data[5], b.data[6], b.data[7]);
    multicomplex3 bv = multicomplex3_make(b.data[8], b.data[9], b.data[10], b.data[11], b.data[12], b.data[13], b.data[14], b.data[15]);
    multicomplex3 p, q;
    multicomplex4 result;
    int i;
    if (k.data[3] == 0) {
        return multicomplex4_exp(k, multicomplex4_mul(k, multicomplex4_log(k, a), b));
    }
    if (k.data[3] > 0) {
        double s = sqrt(k.data[3]);
        p = multicomplex3_pow(base, multicomplex3_add(u, vec8_scale(v, s)), multicomplex3_add(bu, vec8_scale(bv, s)));
        q = multicomplex3_pow(base, multicomplex3_sub(u, vec8_scale(v, s)), multicomplex3_sub(bu, vec8_scale(bv, s)));
        for (i = 0; i < 8; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+8] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[3] < 0)) {
        return multicomplex4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
    }
    /* Find a complex direction in the lower algebra. */
    int bit = 0;
    double inner = 0;
    if (k.data[0] < 0) {
        bit = 1;
        inner = k.data[0];
    }
    else if (k.data[1] < 0) {
        bit = 2;
        inner = k.data[1];
    }
    else if (k.data[2] < 0) {
        bit = 4;
        inner = k.data[2];
    }
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec4 swapped = {{k.data[3], k.data[1], k.data[2], k.data[0]}};
        return multicomplex4_swap_outer(multicomplex4_pow(swapped, multicomplex4_swap_outer(a), multicomplex4_swap_outer(b)));
    }
    double outer_root = sqrt(-k.data[3]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex3 jv;
    multicomplex3 jbv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 8; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
        jbv.data[i] = bv.data[i^bit]*(i & bit ? odd : even);
    }
    p = multicomplex3_pow(base, multicomplex3_add(u, jv), multicomplex3_add(bu, jbv));
    q = multicomplex3_pow(base, multicomplex3_sub(u, jv), multicomplex3_sub(bu, jbv));
    for (i = 0; i < 8; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+8] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex4
multicomplex4_square(vec4 k, multicomplex4 a)
{
    return multicomplex4_mul(k, a, a);
}

static inline multicomplex4
multicomplex4_radd(double a, multicomplex4 b)
{
    return multicomplex4_add(multicomplex4_from_re(a), b);
}

static inline multicomplex4
multicomplex4_ladd(multicomplex4 a, double b)
{
    return multicomplex4_add(a, multicomplex4_from_re(b));
}

static inline multicomplex4
multicomplex4_rmul(vec4 k, double a, multicomplex4 b)
{
    return multicomplex4_mul(k, multicomplex4_from_re(a), b);
}

static inline multicomplex4
multicomplex4_lmul(vec4 k, multicomplex4 a, double b)
{
    return multicomplex4_mul(k, a, multicomplex4_from_re(b));
}

static inline multicomplex4
multicomplex4_rdiv(vec4 k, double a, multicomplex4 b)
{
    return multicomplex4_div(k, multicomplex4_from_re(a), b);
}

static inline multicomplex4
multicomplex4_ldiv(vec4 k, multicomplex4 a, double b)
{
    return multicomplex4_div(k, a, multicomplex4_from_re(b));
}

static inline multicomplex4
multicomplex4_rpow(vec4 k, double a, multicomplex4 b)
{
    return multicomplex4_pow(k, multicomplex4_from_re(a), b);
}

static inline multicomplex4
multicomplex4_lpow(vec4 k, multicomplex4 a, double b)
{
    return multicomplex4_pow(k, a, multicomplex4_from_re(b));
}

static inline multicomplex4
multicomplex4_expm1(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_expm1(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_expm1(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_expm1(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_expm1, multicomplex3_expm1_deriv1);
}

static inline multicomplex4
multicomplex4_log1p(vec4 k, multicomplex4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return complex4_log1p(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1) {
        return splitcomplex4_log1p(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0) {
        return dual4_log1p(a);
    }
    return multicomplex4_apply(k, a, multicomplex3_log1p, multicomplex3_log1p_deriv1);
}

static inline multicomplex4
multicomplex4_mix(vec4 k, multicomplex4 a, multicomplex4 b, double t)
{
    multicomplex4 delta = multicomplex4_mul(k, multicomplex4_inv(k, a), b);
    return multicomplex4_mul(k, a, multicomplex4_lpow(k, delta, t));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
multicomplex4_lmat(vec4 k, multicomplex4 a, double result[ALGPY_MULTICOMPLEX4_N*ALGPY_MULTICOMPLEX4_N])
{
    multicomplex4 coeffs[ALGPY_MULTICOMPLEX4_N * ALGPY_MULTICOMPLEX4_N];
    int row, col, i;
    multicomplex4_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_MULTICOMPLEX4_N; row++) {
        for (col = 0; col < ALGPY_MULTICOMPLEX4_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_MULTICOMPLEX4_N + col].data[row];
            }
            result[row * ALGPY_MULTICOMPLEX4_N + col] = value;
        }
    }
}

static inline void
multicomplex4_rmat(vec4 k, multicomplex4 a, double result[ALGPY_MULTICOMPLEX4_N*ALGPY_MULTICOMPLEX4_N])
{
    multicomplex4 coeffs[ALGPY_MULTICOMPLEX4_N * ALGPY_MULTICOMPLEX4_N];
    int row, col, i;
    multicomplex4_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_MULTICOMPLEX4_N; row++) {
        for (col = 0; col < ALGPY_MULTICOMPLEX4_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_MULTICOMPLEX4_N + i].data[col];
            }
            result[row * ALGPY_MULTICOMPLEX4_N + col] = value;
        }
    }
}


#endif /* ALGPY_MULTICOMPLEX4_H */
