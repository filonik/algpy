/*
 * Core arithmetic of the multicomplex numbers
 * x + x0*e0 + x1*e1 + x01*e0*e1 + x2*e2 + x02*e0*e2 + x12*e1*e2 + x012*e0*e1*e2 with e0*e0 = k[0], e1*e1 = k[1], e2*e2 = k[2].
 *
 * Pure C with no dependency on Python or NumPy. The parameters k are passed
 * explicitly to every operation that depends on them.
 */
#ifndef ALGPY_MULTICOMPLEX3_H
#define ALGPY_MULTICOMPLEX3_H

#include <stdint.h>
#include "../../../fields/complex.h"

#define ALGPY_MULTICOMPLEX3_N 8

#include <math.h>

#include "../../multicomplex/complex3.h"
#include "../../multicomplex/splitcomplex3.h"
#include "../../multicomplex/dual3.h"

#include "multicomplex2.h"
#include "../../../vectorspaces/vec3.h"
#include "../../../vectorspaces/vec8.h"

typedef vec8 multicomplex3;

static inline multicomplex3
multicomplex3_make(double x, double x0, double x1, double x01, double x2, double x02, double x12, double x012)
{
    multicomplex3 r = {{x, x0, x1, x01, x2, x02, x12, x012}};
    return r;
}

static inline multicomplex3
multicomplex3_from_re(double a)
{
    return multicomplex3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
multicomplex3_to_re(multicomplex3 a)
{
    return a.data[0];
}

static inline multicomplex3
multicomplex3_one(void)
{
    return multicomplex3_from_re(1);
}

static inline multicomplex3
multicomplex3_rsub(double a, multicomplex3 b) {
  return vec8_sub(multicomplex3_from_re(a), b);
}

static inline multicomplex3
multicomplex3_lsub(multicomplex3 a, double b) {
  return vec8_sub(a, multicomplex3_from_re(b));
}

static inline multicomplex3
multicomplex3_mul(vec3 h, multicomplex3 a, multicomplex3 b)
{
    const double *k = h.data, *x = a.data, *y = b.data;
    return multicomplex3_make(
      x[0]*y[0] + k[0]*x[1]*y[1] + k[1]*x[2]*y[2] + k[0]*k[1]*x[3]*y[3] + k[2]*x[4]*y[4] + k[0]*k[2]*x[5]*y[5] + k[1]*k[2]*x[6]*y[6] + k[0]*k[1]*k[2]*x[7]*y[7],
      x[0]*y[1] + k[1]*k[2]*x[6]*y[7] + k[1]*k[2]*x[7]*y[6] + k[1]*x[2]*y[3] + k[1]*x[3]*y[2] + k[2]*x[4]*y[5] + k[2]*x[5]*y[4] + x[1]*y[0],
      x[0]*y[2] + k[0]*k[2]*x[5]*y[7] + k[0]*k[2]*x[7]*y[5] + k[0]*x[1]*y[3] + k[0]*x[3]*y[1] + k[2]*x[4]*y[6] + k[2]*x[6]*y[4] + x[2]*y[0],
      x[0]*y[3] + k[2]*x[4]*y[7] + k[2]*x[5]*y[6] + k[2]*x[6]*y[5] + k[2]*x[7]*y[4] + x[1]*y[2] + x[2]*y[1] + x[3]*y[0],
      x[0]*y[4] + k[0]*k[1]*x[3]*y[7] + k[0]*k[1]*x[7]*y[3] + k[0]*x[1]*y[5] + k[0]*x[5]*y[1] + k[1]*x[2]*y[6] + k[1]*x[6]*y[2] + x[4]*y[0],
      x[0]*y[5] + k[1]*x[2]*y[7] + k[1]*x[3]*y[6] + k[1]*x[6]*y[3] + k[1]*x[7]*y[2] + x[1]*y[4] + x[4]*y[1] + x[5]*y[0],
      x[0]*y[6] + k[0]*x[1]*y[7] + k[0]*x[3]*y[5] + k[0]*x[5]*y[3] + k[0]*x[7]*y[1] + x[2]*y[4] + x[4]*y[2] + x[6]*y[0],
      x[0]*y[7] + x[1]*y[6] + x[2]*y[5] + x[3]*y[4] + x[4]*y[3] + x[5]*y[2] + x[6]*y[1] + x[7]*y[0]);
}

static inline multicomplex3
multicomplex3_conjugate(multicomplex3 a)
{
    return multicomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3],
                              -a.data[4], -a.data[5], -a.data[6], -a.data[7]);
}

/* Norm a*conj(a) = u^2 - k[2]*v^2 with a = u + v*e2. */
static inline multicomplex2
multicomplex3_norm2(vec3 k, multicomplex3 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    multicomplex2 u = multicomplex2_make(a.data[0], a.data[1],
                                         a.data[2], a.data[3]);
    multicomplex2 v = multicomplex2_make(a.data[4], a.data[5],
                                         a.data[6], a.data[7]);
    multicomplex2 uu = multicomplex2_mul(k01, u, u);
    multicomplex2 vv = multicomplex2_mul(k01, v, v);
    return vec4_sub(uu, vec4_scale(vv, k.data[2]));
}

static inline multicomplex1
multicomplex3_norm4(vec3 k, multicomplex3 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    return multicomplex2_norm2(k01, multicomplex3_norm2(k, a));
}

static inline double
multicomplex3_norm8(vec3 k, multicomplex3 a)
{
    vec1 k0 = {{k.data[0]}};
    return multicomplex1_norm2(k0, multicomplex3_norm4(k, a));
}

static inline complex
multicomplex3_norm(vec3 k, multicomplex3 a)
{
    return complex_root_real(multicomplex3_norm8(k, a), 8);
}

static inline multicomplex3
multicomplex3_inv(vec3 k, multicomplex3 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    multicomplex2 ninv = multicomplex2_inv(k01, multicomplex3_norm2(k, a));
    multicomplex2 u = multicomplex2_make(a.data[0], a.data[1],
                                         a.data[2], a.data[3]);
    multicomplex2 v = multicomplex2_make(a.data[4], a.data[5],
                                         a.data[6], a.data[7]);
    multicomplex2 ru = multicomplex2_mul(k01, u, ninv);
    multicomplex2 rv = multicomplex2_mul(k01, vec4_neg(v), ninv);
    return multicomplex3_make(ru.data[0], ru.data[1], ru.data[2], ru.data[3],
                              rv.data[0], rv.data[1], rv.data[2], rv.data[3]);
}

static inline multicomplex3
multicomplex3_div(vec3 k, multicomplex3 a, multicomplex3 b)
{
    if (vec8_is_real(b)) {
        return vec8_div_scalar(a, b.data[0]);
    }
    return multicomplex3_mul(k, a, multicomplex3_inv(k, b));
}

static inline void
multicomplex3_basis(multicomplex3 result[ALGPY_MULTICOMPLEX3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_MULTICOMPLEX3_N; i++) {
        double components[ALGPY_MULTICOMPLEX3_N] = {0};
        components[i] = 1;
        for (j = 0; j < ALGPY_MULTICOMPLEX3_N; j++) {
            result[i].data[j] = components[j];
        }
    }
}

static inline void
multicomplex3_coeffs(
    vec3 k,
    multicomplex3 result[ALGPY_MULTICOMPLEX3_N * ALGPY_MULTICOMPLEX3_N])
{
    multicomplex3 basis[ALGPY_MULTICOMPLEX3_N];
    int i, j;
    multicomplex3_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX3_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX3_N; j++) {
            result[i * ALGPY_MULTICOMPLEX3_N + j] =
                multicomplex3_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
multicomplex3_minpoly(vec3 k, multicomplex3 a, double result[9])
{
    const int n = ALGPY_MULTICOMPLEX3_N;
    multicomplex3 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * multicomplex3_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            multicomplex3 adj = multicomplex3_lsub(uk, ck);
            uk = multicomplex3_mul(k, a, adj);
        }
    }
}

/*
 * Eigenvalues sum of x_S * prod_{i in S} (+/- sqrt(k_i)) over all subsets S,
 * with the sign of sqrt(k_i) given by the bits of the eigenvalue index,
 * as complex numbers.
 */
static inline void
multicomplex3_eivals(vec3 k, multicomplex3 a, complex result[8])
{
    complex s[3];
    complex n[8];
    int i, j, b;

    for (b = 0; b < 3; b++) {
        s[b] = complex_sqrt_real(k.data[b]);
    }
    for (i = 0; i < 8; i++) {
        complex t = complex_make(a.data[i], 0);
        for (b = 0; b < 3; b++) {
            if (i & (1 << b)) {
                t = complex_mul(t, s[b]);
            }
        }
        n[i] = t;
    }
    for (j = 0; j < 8; j++) {
        double re = 0, im = 0;
        for (i = 0; i < 8; i++) {
            int sign = 1;
            for (b = 0; b < 3; b++) {
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
static inline multicomplex3
multicomplex3_zero(void)
{
    return vec8_zero();
}

static inline int
multicomplex3_is_real(multicomplex3 a)
{
    return vec8_is_real(a);
}

static inline multicomplex3
multicomplex3_pos(multicomplex3 a)
{
    return vec8_pos(a);
}

static inline multicomplex3
multicomplex3_neg(multicomplex3 a)
{
    return vec8_neg(a);
}

static inline multicomplex3
multicomplex3_add(multicomplex3 a, multicomplex3 b)
{
    return vec8_add(a, b);
}

static inline multicomplex3
multicomplex3_sub(multicomplex3 a, multicomplex3 b)
{
    return vec8_sub(a, b);
}

static inline int
multicomplex3_equal(multicomplex3 a, multicomplex3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
multicomplex3_not_equal(multicomplex3 a, multicomplex3 b)
{
    return !multicomplex3_equal(a, b);
}

static inline multicomplex3
multicomplex3_re(multicomplex3 a)
{
    return multicomplex3_from_re(a.data[0]);
}

static inline multicomplex3
multicomplex3_im(multicomplex3 a)
{
    return multicomplex3_sub(a, multicomplex3_re(a));
}

/* Exchange the first and outer generators, which commute. */
static inline multicomplex3
multicomplex3_swap_outer(multicomplex3 a)
{
    return multicomplex3_make(a.data[0], a.data[4], a.data[2], a.data[6], a.data[1], a.data[5], a.data[3], a.data[7]);
}

static inline multicomplex3
multicomplex3_apply(vec3 k, multicomplex3 a, multicomplex2 (*f)(vec2, multicomplex2), multicomplex2 (*deriv1)(vec2, multicomplex2))
{
    vec2 base = {{k.data[0], k.data[1]}};
    multicomplex2 u = multicomplex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    multicomplex2 v = multicomplex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    multicomplex2 p, q;
    multicomplex3 result;
    int i;
    if (k.data[2] == 0) {
        p = f(base, u);
        q = multicomplex2_mul(base, v, deriv1(base, u));
        for (i = 0; i < 4; i++) {
            result.data[i] = p.data[i];
            result.data[i+4] = q.data[i];
        }
        return result;
    }
    if (k.data[2] > 0) {
        double s = sqrt(k.data[2]);
        p = f(base, multicomplex2_add(u, vec4_scale(v, s)));
        q = f(base, multicomplex2_sub(u, vec4_scale(v, s)));
        for (i = 0; i < 4; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+4] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[2] < 0)) {
        return multicomplex3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
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
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec3 swapped = {{k.data[2], k.data[1], k.data[0]}};
        return multicomplex3_swap_outer(multicomplex3_apply(swapped, multicomplex3_swap_outer(a), f, deriv1));
    }
    double outer_root = sqrt(-k.data[2]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex2 jv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 4; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
    }
    p = f(base, multicomplex2_add(u, jv));
    q = f(base, multicomplex2_sub(u, jv));
    for (i = 0; i < 4; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+4] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex3
multicomplex3_exp(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_exp(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_exp(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_exp(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_exp, multicomplex2_exp_deriv1);
}

static inline multicomplex3
multicomplex3_log(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_log(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_log(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_log(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_log, multicomplex2_log_deriv1);
}

static inline multicomplex3
multicomplex3_sqrt(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_sqrt(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_sqrt(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_sqrt(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_sqrt, multicomplex2_sqrt_deriv1);
}

static inline multicomplex3
multicomplex3_cos(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_cos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_cos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_cos(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_cos, multicomplex2_cos_deriv1);
}

static inline multicomplex3
multicomplex3_sin(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_sin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_sin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_sin(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_sin, multicomplex2_sin_deriv1);
}

static inline multicomplex3
multicomplex3_tan(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_tan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_tan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_tan(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_tan, multicomplex2_tan_deriv1);
}

static inline multicomplex3
multicomplex3_acos(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_acos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_acos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_acos(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_acos, multicomplex2_acos_deriv1);
}

static inline multicomplex3
multicomplex3_asin(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_asin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_asin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_asin(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_asin, multicomplex2_asin_deriv1);
}

static inline multicomplex3
multicomplex3_atan(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_atan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_atan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_atan(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_atan, multicomplex2_atan_deriv1);
}

static inline multicomplex3
multicomplex3_cosh(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_cosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_cosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_cosh(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_cosh, multicomplex2_cosh_deriv1);
}

static inline multicomplex3
multicomplex3_sinh(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_sinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_sinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_sinh(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_sinh, multicomplex2_sinh_deriv1);
}

static inline multicomplex3
multicomplex3_tanh(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_tanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_tanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_tanh(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_tanh, multicomplex2_tanh_deriv1);
}

static inline multicomplex3
multicomplex3_acosh(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_acosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_acosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_acosh(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_acosh, multicomplex2_acosh_deriv1);
}

static inline multicomplex3
multicomplex3_asinh(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_asinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_asinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_asinh(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_asinh, multicomplex2_asinh_deriv1);
}

static inline multicomplex3
multicomplex3_atanh(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_atanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_atanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_atanh(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_atanh, multicomplex2_atanh_deriv1);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline multicomplex3
multicomplex3_ipow(vec3 k, multicomplex3 a, int64_t b)
{
    multicomplex3 result = multicomplex3_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return multicomplex3_mul(k, a, a);
    }
    if (b == -1) {
        return multicomplex3_inv(k, a);
    }
    if (b < 0) {
        a = multicomplex3_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = multicomplex3_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = multicomplex3_mul(k, a, a);
        }
    }
    return result;
}

static inline multicomplex3
multicomplex3_rpow(vec3 k, multicomplex3 a, multicomplex3 b)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_rpow(a, b);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_rpow(a, b);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_rpow(a, b);
    }
    vec2 base = {{k.data[0], k.data[1]}};
    multicomplex2 u = multicomplex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    multicomplex2 v = multicomplex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    multicomplex2 bu = multicomplex2_make(b.data[0], b.data[1], b.data[2], b.data[3]);
    multicomplex2 bv = multicomplex2_make(b.data[4], b.data[5], b.data[6], b.data[7]);
    multicomplex2 p, q;
    multicomplex3 result;
    int i;
    if (k.data[2] == 0) {
        return multicomplex3_exp(k, multicomplex3_mul(k, multicomplex3_log(k, a), b));
    }
    if (k.data[2] > 0) {
        double s = sqrt(k.data[2]);
        p = multicomplex2_rpow(base, multicomplex2_add(u, vec4_scale(v, s)), multicomplex2_add(bu, vec4_scale(bv, s)));
        q = multicomplex2_rpow(base, multicomplex2_sub(u, vec4_scale(v, s)), multicomplex2_sub(bu, vec4_scale(bv, s)));
        for (i = 0; i < 4; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+4] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[2] < 0)) {
        return multicomplex3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
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
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec3 swapped = {{k.data[2], k.data[1], k.data[0]}};
        return multicomplex3_swap_outer(multicomplex3_rpow(swapped, multicomplex3_swap_outer(a), multicomplex3_swap_outer(b)));
    }
    double outer_root = sqrt(-k.data[2]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex2 jv;
    multicomplex2 jbv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 4; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
        jbv.data[i] = bv.data[i^bit]*(i & bit ? odd : even);
    }
    p = multicomplex2_rpow(base, multicomplex2_add(u, jv), multicomplex2_add(bu, jbv));
    q = multicomplex2_rpow(base, multicomplex2_sub(u, jv), multicomplex2_sub(bu, jbv));
    for (i = 0; i < 4; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+4] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex3
multicomplex3_pow(vec3 k, multicomplex3 a, multicomplex3 b)
{
    return multicomplex3_rpow(k, a, b);
}

static inline multicomplex3
multicomplex3_lpow(vec3 k, multicomplex3 a, multicomplex3 b)
{
    return multicomplex3_rpow(k, a, b);
}

static inline multicomplex3
multicomplex3_square(vec3 k, multicomplex3 a)
{
    return multicomplex3_mul(k, a, a);
}

static inline multicomplex3
multicomplex3_radd(double a, multicomplex3 b)
{
    return multicomplex3_add(multicomplex3_from_re(a), b);
}

static inline multicomplex3
multicomplex3_ladd(multicomplex3 a, double b)
{
    return multicomplex3_add(a, multicomplex3_from_re(b));
}

static inline multicomplex3
multicomplex3_rmul(vec3 k, double a, multicomplex3 b)
{
    return multicomplex3_mul(k, multicomplex3_from_re(a), b);
}

static inline multicomplex3
multicomplex3_lmul(vec3 k, multicomplex3 a, double b)
{
    return multicomplex3_mul(k, a, multicomplex3_from_re(b));
}

static inline multicomplex3
multicomplex3_rdiv(vec3 k, double a, multicomplex3 b)
{
    return multicomplex3_div(k, multicomplex3_from_re(a), b);
}

static inline multicomplex3
multicomplex3_ldiv(vec3 k, multicomplex3 a, double b)
{
    return multicomplex3_div(k, a, multicomplex3_from_re(b));
}

static inline multicomplex3
multicomplex3_rpows(vec3 k, double a, multicomplex3 b)
{
    return multicomplex3_rpow(k, multicomplex3_from_re(a), b);
}

static inline multicomplex3
multicomplex3_lpows(vec3 k, multicomplex3 a, double b)
{
    return multicomplex3_rpow(k, a, multicomplex3_from_re(b));
}

static inline multicomplex3
multicomplex3_expm1(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_expm1(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_expm1(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_expm1(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_expm1, multicomplex2_expm1_deriv1);
}

static inline multicomplex3
multicomplex3_log1p(vec3 k, multicomplex3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return complex3_log1p(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1) {
        return splitcomplex3_log1p(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0) {
        return dual3_log1p(a);
    }
    return multicomplex3_apply(k, a, multicomplex2_log1p, multicomplex2_log1p_deriv1);
}

static inline multicomplex3
multicomplex3_mix(vec3 k, multicomplex3 a, multicomplex3 b, double t)
{
    multicomplex3 delta = multicomplex3_mul(k, multicomplex3_inv(k, a), b);
    return multicomplex3_mul(k, a, multicomplex3_lpows(k, delta, t));
}

/* Derivatives through order 1, for nilpotent extensions at higher levels. */
static inline multicomplex3
multicomplex3_exp_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_exp_deriv1, multicomplex2_exp_deriv2);
}

static inline multicomplex3
multicomplex3_expm1_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_expm1_deriv1, multicomplex2_expm1_deriv2);
}

static inline multicomplex3
multicomplex3_log_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_log_deriv1, multicomplex2_log_deriv2);
}

static inline multicomplex3
multicomplex3_log1p_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_log1p_deriv1, multicomplex2_log1p_deriv2);
}

static inline multicomplex3
multicomplex3_inv_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_inv_deriv1, multicomplex2_inv_deriv2);
}

static inline multicomplex3
multicomplex3_square_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_square_deriv1, multicomplex2_square_deriv2);
}

static inline multicomplex3
multicomplex3_sqrt_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_sqrt_deriv1, multicomplex2_sqrt_deriv2);
}

static inline multicomplex3
multicomplex3_sin_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_sin_deriv1, multicomplex2_sin_deriv2);
}

static inline multicomplex3
multicomplex3_cos_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_cos_deriv1, multicomplex2_cos_deriv2);
}

static inline multicomplex3
multicomplex3_tan_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_tan_deriv1, multicomplex2_tan_deriv2);
}

static inline multicomplex3
multicomplex3_asin_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_asin_deriv1, multicomplex2_asin_deriv2);
}

static inline multicomplex3
multicomplex3_acos_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_acos_deriv1, multicomplex2_acos_deriv2);
}

static inline multicomplex3
multicomplex3_atan_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_atan_deriv1, multicomplex2_atan_deriv2);
}

static inline multicomplex3
multicomplex3_sinh_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_sinh_deriv1, multicomplex2_sinh_deriv2);
}

static inline multicomplex3
multicomplex3_cosh_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_cosh_deriv1, multicomplex2_cosh_deriv2);
}

static inline multicomplex3
multicomplex3_tanh_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_tanh_deriv1, multicomplex2_tanh_deriv2);
}

static inline multicomplex3
multicomplex3_asinh_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_asinh_deriv1, multicomplex2_asinh_deriv2);
}

static inline multicomplex3
multicomplex3_acosh_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_acosh_deriv1, multicomplex2_acosh_deriv2);
}

static inline multicomplex3
multicomplex3_atanh_deriv1(vec3 k, multicomplex3 a)
{
    return multicomplex3_apply(k, a, multicomplex2_atanh_deriv1, multicomplex2_atanh_deriv2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
multicomplex3_lmat(vec3 k, multicomplex3 a, double result[ALGPY_MULTICOMPLEX3_N*ALGPY_MULTICOMPLEX3_N])
{
    multicomplex3 basis = {{0}};
    multicomplex3 product;
    int row, col;
    for (col = 0; col < ALGPY_MULTICOMPLEX3_N; col++) {
        basis.data[col] = 1;
        product = multicomplex3_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_MULTICOMPLEX3_N; row++) {
            result[row * ALGPY_MULTICOMPLEX3_N + col] = product.data[row];
        }
    }
}

static inline void
multicomplex3_rmat(vec3 k, multicomplex3 a, double result[ALGPY_MULTICOMPLEX3_N*ALGPY_MULTICOMPLEX3_N])
{
    multicomplex3 basis = {{0}};
    multicomplex3 product;
    int row, col;
    for (row = 0; row < ALGPY_MULTICOMPLEX3_N; row++) {
        basis.data[row] = 1;
        product = multicomplex3_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_MULTICOMPLEX3_N; col++) {
            result[row * ALGPY_MULTICOMPLEX3_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_MULTICOMPLEX3_H */
