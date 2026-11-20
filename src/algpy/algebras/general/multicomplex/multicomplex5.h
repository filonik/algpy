/*
 * Core arithmetic of the multicomplex numbers
 * sum_S x_S * e_S over subsets S of {0, 1, 2, 3, 4} (component index = bitmask of S)
 * with e0*e0 = k[0], e1*e1 = k[1], e2*e2 = k[2], e3*e3 = k[3], e4*e4 = k[4].
 *
 * Pure C with no dependency on Python or NumPy. The parameters k are passed
 * explicitly to every operation that depends on them.
 */
#ifndef ALGPY_MULTICOMPLEX5_H
#define ALGPY_MULTICOMPLEX5_H

#include <string.h>
#include <stdint.h>
#include "../../../fields/complex.h"

#define ALGPY_MULTICOMPLEX5_N 32

#include <math.h>

#include "../../multicomplex/complex5.h"
#include "../../multicomplex/splitcomplex5.h"
#include "../../multicomplex/dual5.h"

#include "multicomplex4.h"
#include "../../../vectorspaces/vec5.h"
#include "../../../vectorspaces/vec32.h"

typedef vec32 multicomplex5;

static inline multicomplex5
multicomplex5_make(const double coeffs[ALGPY_MULTICOMPLEX5_N])
{
    multicomplex5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline multicomplex5
multicomplex5_from_re(double a)
{
    multicomplex5 result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
multicomplex5_to_re(multicomplex5 a)
{
    return a.data[0];
}

static inline multicomplex5
multicomplex5_one(void)
{
    return multicomplex5_from_re(1);
}

static inline multicomplex5
multicomplex5_rsub(double a, multicomplex5 b)
{
    return vec32_sub(multicomplex5_from_re(a), b);
}

static inline multicomplex5
multicomplex5_lsub(multicomplex5 a, double b)
{
    return vec32_sub(a, multicomplex5_from_re(b));
}

/* e_i * e_j = (product of k[b] over the bits b shared by i and j) * e_(i^j). */
static inline multicomplex5
multicomplex5_mul(vec5 h, multicomplex5 a, multicomplex5 b)
{
    const double *k = h.data;
    multicomplex5 r;
    int i, j, c;
    for (i = 0; i < ALGPY_MULTICOMPLEX5_N; i++) {
        r.data[i] = 0;
    }
    for (i = 0; i < ALGPY_MULTICOMPLEX5_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX5_N; j++) {
            double t = a.data[i] * b.data[j];
            for (c = 0; c < 5; c++) {
                if ((i & j) & (1 << c)) {
                    t *= k[c];
                }
            }
            r.data[i ^ j] += t;
        }
    }
    return r;
}

static inline multicomplex5
multicomplex5_conjugate(multicomplex5 a)
{
    multicomplex5 r = a;
    int i;
    for (i = 16; i < ALGPY_MULTICOMPLEX5_N; i++) {
        r.data[i] = -r.data[i];
    }
    return r;
}

static inline multicomplex4
multicomplex5_lower(multicomplex5 a)
{
    return multicomplex4_make(a.data);
}

static inline multicomplex4
multicomplex5_upper(multicomplex5 a)
{
    return multicomplex4_make(a.data + 16);
}

static inline multicomplex5
multicomplex5_join(multicomplex4 u, multicomplex4 v)
{
    multicomplex5 result;
    memcpy(result.data, u.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_MULTICOMPLEX5_N/2, v.data, sizeof(result.data)/2);
    return result;
}

/* Norm a*conj(a) = u^2 - k[4]*v^2 with a = u + v*e4. */
static inline multicomplex4
multicomplex5_norm2(vec5 k, multicomplex5 a)
{
    vec4 base = {{k.data[0], k.data[1], k.data[2], k.data[3]}};
    multicomplex4 u = multicomplex5_lower(a);
    multicomplex4 v = multicomplex5_upper(a);
    return vec16_sub(multicomplex4_mul(base, u, u),
                     vec16_scale(multicomplex4_mul(base, v, v), k.data[4]));
}

static inline multicomplex3
multicomplex5_norm4(vec5 k, multicomplex5 a)
{
    vec4 base = {{k.data[0], k.data[1], k.data[2], k.data[3]}};
    return multicomplex4_norm2(base, multicomplex5_norm2(k, a));
}

static inline multicomplex2
multicomplex5_norm8(vec5 k, multicomplex5 a)
{
    vec3 base = {{k.data[0], k.data[1], k.data[2]}};
    return multicomplex3_norm2(base, multicomplex5_norm4(k, a));
}

static inline multicomplex1
multicomplex5_norm16(vec5 k, multicomplex5 a)
{
    vec2 base = {{k.data[0], k.data[1]}};
    return multicomplex2_norm2(base, multicomplex5_norm8(k, a));
}

static inline double
multicomplex5_norm32(vec5 k, multicomplex5 a)
{
    vec1 base = {{k.data[0]}};
    return multicomplex1_norm2(base, multicomplex5_norm16(k, a));
}

static inline complex
multicomplex5_norm(vec5 k, multicomplex5 a)
{
    return complex_root_real(multicomplex5_norm32(k, a), 32);
}

static inline multicomplex5
multicomplex5_inv(vec5 k, multicomplex5 a)
{
    vec4 k0123 = {{k.data[0], k.data[1], k.data[2], k.data[3]}};
    multicomplex4 ninv = multicomplex4_inv(k0123, multicomplex5_norm2(k, a));
    multicomplex4 u = multicomplex5_lower(a);
    multicomplex4 v = multicomplex5_upper(a);
    multicomplex4 ru = multicomplex4_mul(k0123, u, ninv);
    multicomplex4 rv = multicomplex4_mul(k0123, vec16_neg(v), ninv);
    return multicomplex5_join(ru, rv);
}

static inline multicomplex5
multicomplex5_div(vec5 k, multicomplex5 a, multicomplex5 b)
{
    if (vec32_is_real(b)) {
        return vec32_div_scalar(a, b.data[0]);
    }
    return multicomplex5_mul(k, a, multicomplex5_inv(k, b));
}

static inline void
multicomplex5_basis(multicomplex5 result[ALGPY_MULTICOMPLEX5_N])
{
    int i, j;
    for (i = 0; i < ALGPY_MULTICOMPLEX5_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX5_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
multicomplex5_gens(multicomplex5 result[5])
{
    result[0] = multicomplex5_make((const double[ALGPY_MULTICOMPLEX5_N]){[1] = 1});
    result[1] = multicomplex5_make((const double[ALGPY_MULTICOMPLEX5_N]){[2] = 1});
    result[2] = multicomplex5_make((const double[ALGPY_MULTICOMPLEX5_N]){[4] = 1});
    result[3] = multicomplex5_make((const double[ALGPY_MULTICOMPLEX5_N]){[8] = 1});
    result[4] = multicomplex5_make((const double[ALGPY_MULTICOMPLEX5_N]){[16] = 1});
}

static inline void
multicomplex5_coeffs(
    vec5 k,
    multicomplex5 result[ALGPY_MULTICOMPLEX5_N * ALGPY_MULTICOMPLEX5_N])
{
    multicomplex5 basis[ALGPY_MULTICOMPLEX5_N];
    int i, j;
    multicomplex5_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX5_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX5_N; j++) {
            result[i * ALGPY_MULTICOMPLEX5_N + j] =
                multicomplex5_mul(k, basis[i], basis[j]);
        }
    }
}

/*
 * Eigenvalues sum of x_S * prod_{i in S} (+/- sqrt(k_i)) over all subsets S,
 * with the sign of sqrt(k_i) given by the bits of the eigenvalue index,
 * as complex numbers.
 */
static inline void
multicomplex5_eivals(vec5 k, multicomplex5 a, complex result[32])
{
    complex s[5];
    complex n[32];
    int i, j, b;

    for (b = 0; b < 5; b++) {
        s[b] = complex_sqrt_real(k.data[b]);
    }
    for (i = 0; i < 32; i++) {
        complex t = complex_make((const double[]){a.data[i], 0});
        for (b = 0; b < 5; b++) {
            if (i & (1 << b)) {
                t = complex_mul(t, s[b]);
            }
        }
        n[i] = t;
    }
    for (j = 0; j < 32; j++) {
        double re = 0, im = 0;
        for (i = 0; i < 32; i++) {
            int sign = 1;
            for (b = 0; b < 5; b++) {
                if ((i & (1 << b)) && !(j & (1 << b))) {
                    sign = -sign;
                }
            }
            re += sign * n[i].data[0];
            im += sign * n[i].data[1];
        }
        result[j] = complex_make((const double[]){re, im});
    }
}

/* Form the generic degree-32 polynomial from the explicit spectrum.
 * Pair conjugate roots into real quadratics to keep all coefficients real
 * and avoid cancellation in the higher-degree trace recurrence.
 * Zero parameters retain the repeated roots of the generic polynomial. */
static inline void
multicomplex5_minpoly(vec5 k, multicomplex5 a, double result[33])
{
    complex roots[32];
    int bit = 0, degree = 0;
    int i, j;
    multicomplex5_eivals(k, a, roots);
    for (i = 0; i < 5; i++) {
        if (k.data[i] < 0) {
            bit = 1 << i;
            break;
        }
    }
    for (i = 0; i <= 32; i++) result[i] = 0;
    result[0] = 1;
    for (i = 0; i < 32; i++) {
        if (bit) {
            if (i & bit) continue;
            double linear = -2*roots[i].data[0];
            double constant = complex_norm2(roots[i]);
            for (j = degree+2; j >= 2; j--) {
                result[j] = result[j-2] + linear*result[j-1] + constant*result[j];
            }
            result[1] = linear*result[0] + constant*result[1];
            result[0] *= constant;
            degree += 2;
        }
        else {
            double root = roots[i].data[0];
            for (j = degree+1; j > 0; j--) {
                result[j] = result[j-1] - root*result[j];
            }
            result[0] *= -root;
            degree++;
        }
    }
}

/* Common algebra API. */
static inline multicomplex5
multicomplex5_zero(void)
{
    return vec32_zero();
}

static inline int
multicomplex5_is_real(multicomplex5 a)
{
    return vec32_is_real(a);
}

static inline multicomplex5
multicomplex5_pos(multicomplex5 a)
{
    return vec32_pos(a);
}

static inline multicomplex5
multicomplex5_neg(multicomplex5 a)
{
    return vec32_neg(a);
}

static inline multicomplex5
multicomplex5_add(multicomplex5 a, multicomplex5 b)
{
    return vec32_add(a, b);
}

static inline multicomplex5
multicomplex5_sub(multicomplex5 a, multicomplex5 b)
{
    return vec32_sub(a, b);
}

static inline int
multicomplex5_equal(multicomplex5 a, multicomplex5 b)
{
    return vec32_all_eq(a, b);
}

static inline int
multicomplex5_not_equal(multicomplex5 a, multicomplex5 b)
{
    return !multicomplex5_equal(a, b);
}

static inline multicomplex5
multicomplex5_re(multicomplex5 a)
{
    return multicomplex5_from_re(a.data[0]);
}

static inline multicomplex5
multicomplex5_im(multicomplex5 a)
{
    return multicomplex5_sub(a, multicomplex5_re(a));
}

/* Exchange the first and outer generators, which commute. */
static inline multicomplex5
multicomplex5_swap_outer(multicomplex5 a)
{
    return multicomplex5_make((const double[]){a.data[0], a.data[16], a.data[2], a.data[18], a.data[4], a.data[20], a.data[6], a.data[22], a.data[8], a.data[24], a.data[10], a.data[26], a.data[12], a.data[28], a.data[14], a.data[30], a.data[1], a.data[17], a.data[3], a.data[19], a.data[5], a.data[21], a.data[7], a.data[23], a.data[9], a.data[25], a.data[11], a.data[27], a.data[13], a.data[29], a.data[15], a.data[31]});
}

static inline multicomplex5
multicomplex5_apply(vec5 k, multicomplex5 a, multicomplex4 (*f)(vec4, multicomplex4), multicomplex4 (*deriv1)(vec4, multicomplex4))
{
    vec4 base = {{k.data[0], k.data[1], k.data[2], k.data[3]}};
    multicomplex4 u = multicomplex5_lower(a);
    multicomplex4 v = multicomplex5_upper(a);
    multicomplex4 p, q;
    multicomplex5 result;
    int i;
    if (k.data[4] == 0) {
        p = f(base, u);
        q = multicomplex4_mul(base, v, deriv1(base, u));
        for (i = 0; i < 16; i++) {
            result.data[i] = p.data[i];
            result.data[i+16] = q.data[i];
        }
        return result;
    }
    if (k.data[4] > 0) {
        double s = sqrt(k.data[4]);
        p = f(base, multicomplex4_add(u, vec16_scale(v, s)));
        q = f(base, multicomplex4_sub(u, vec16_scale(v, s)));
        for (i = 0; i < 16; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+16] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[4] < 0)) {
        return multicomplex5_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
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
    else if (k.data[3] < 0) {
        bit = 8;
        inner = k.data[3];
    }
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec5 swapped = {{k.data[4], k.data[1], k.data[2], k.data[3], k.data[0]}};
        return multicomplex5_swap_outer(multicomplex5_apply(swapped, multicomplex5_swap_outer(a), f, deriv1));
    }
    double outer_root = sqrt(-k.data[4]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex4 jv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 16; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
    }
    p = f(base, multicomplex4_add(u, jv));
    q = f(base, multicomplex4_sub(u, jv));
    for (i = 0; i < 16; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+16] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex5
multicomplex5_exp(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_exp(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_exp(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_exp(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_exp, multicomplex4_exp_deriv1);
}

static inline multicomplex5
multicomplex5_log(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_log(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_log(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_log(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_log, multicomplex4_log_deriv1);
}

static inline multicomplex5
multicomplex5_sqrt(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_sqrt(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_sqrt(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_sqrt(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_sqrt, multicomplex4_sqrt_deriv1);
}

static inline multicomplex5
multicomplex5_cos(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_cos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_cos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_cos(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_cos, multicomplex4_cos_deriv1);
}

static inline multicomplex5
multicomplex5_sin(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_sin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_sin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_sin(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_sin, multicomplex4_sin_deriv1);
}

static inline multicomplex5
multicomplex5_tan(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_tan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_tan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_tan(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_tan, multicomplex4_tan_deriv1);
}

static inline multicomplex5
multicomplex5_acos(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_acos(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_acos(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_acos(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_acos, multicomplex4_acos_deriv1);
}

static inline multicomplex5
multicomplex5_asin(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_asin(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_asin(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_asin(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_asin, multicomplex4_asin_deriv1);
}

static inline multicomplex5
multicomplex5_atan(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_atan(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_atan(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_atan(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_atan, multicomplex4_atan_deriv1);
}

static inline multicomplex5
multicomplex5_cosh(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_cosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_cosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_cosh(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_cosh, multicomplex4_cosh_deriv1);
}

static inline multicomplex5
multicomplex5_sinh(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_sinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_sinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_sinh(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_sinh, multicomplex4_sinh_deriv1);
}

static inline multicomplex5
multicomplex5_tanh(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_tanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_tanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_tanh(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_tanh, multicomplex4_tanh_deriv1);
}

static inline multicomplex5
multicomplex5_acosh(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_acosh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_acosh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_acosh(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_acosh, multicomplex4_acosh_deriv1);
}

static inline multicomplex5
multicomplex5_asinh(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_asinh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_asinh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_asinh(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_asinh, multicomplex4_asinh_deriv1);
}

static inline multicomplex5
multicomplex5_atanh(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_atanh(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_atanh(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_atanh(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_atanh, multicomplex4_atanh_deriv1);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline multicomplex5
multicomplex5_ipow(vec5 k, multicomplex5 a, int64_t b)
{
    multicomplex5 result = multicomplex5_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return multicomplex5_mul(k, a, a);
    }
    if (b == -1) {
        return multicomplex5_inv(k, a);
    }
    if (b < 0) {
        a = multicomplex5_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = multicomplex5_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = multicomplex5_mul(k, a, a);
        }
    }
    return result;
}

static inline multicomplex5
multicomplex5_rpow(vec5 k, multicomplex5 a, multicomplex5 b)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_rpow(a, b);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_rpow(a, b);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_rpow(a, b);
    }
    vec4 base = {{k.data[0], k.data[1], k.data[2], k.data[3]}};
    multicomplex4 u = multicomplex5_lower(a);
    multicomplex4 v = multicomplex5_upper(a);
    multicomplex4 bu = multicomplex5_lower(b);
    multicomplex4 bv = multicomplex5_upper(b);
    multicomplex4 p, q;
    multicomplex5 result;
    int i;
    if (k.data[4] == 0) {
        return multicomplex5_exp(k, multicomplex5_mul(k, multicomplex5_log(k, a), b));
    }
    if (k.data[4] > 0) {
        double s = sqrt(k.data[4]);
        p = multicomplex4_rpow(base, multicomplex4_add(u, vec16_scale(v, s)), multicomplex4_add(bu, vec16_scale(bv, s)));
        q = multicomplex4_rpow(base, multicomplex4_sub(u, vec16_scale(v, s)), multicomplex4_sub(bu, vec16_scale(bv, s)));
        for (i = 0; i < 16; i++) {
            result.data[i] = (p.data[i]+q.data[i])/2;
            result.data[i+16] = ((p.data[i]-q.data[i])/s)/2;
        }
        return result;
    }
    if (!(k.data[4] < 0)) {
        return multicomplex5_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
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
    else if (k.data[3] < 0) {
        bit = 8;
        inner = k.data[3];
    }
    else {
        /* Bring the first negative generator into the lower algebra. */
        vec5 swapped = {{k.data[4], k.data[1], k.data[2], k.data[3], k.data[0]}};
        return multicomplex5_swap_outer(multicomplex5_rpow(swapped, multicomplex5_swap_outer(a), multicomplex5_swap_outer(b)));
    }
    double outer_root = sqrt(-k.data[4]);
    double inner_root = sqrt(-inner);
    double even = -outer_root*inner_root;
    double odd = outer_root/inner_root;
    multicomplex4 jv;
    multicomplex4 jbv;
    /* j = sqrt(k_outer/k_inner)*e_inner, so j*j = k_outer. */
    for (i = 0; i < 16; i++) {
        jv.data[i] = v.data[i^bit]*(i & bit ? odd : even);
        jbv.data[i] = bv.data[i^bit]*(i & bit ? odd : even);
    }
    p = multicomplex4_rpow(base, multicomplex4_add(u, jv), multicomplex4_add(bu, jbv));
    q = multicomplex4_rpow(base, multicomplex4_sub(u, jv), multicomplex4_sub(bu, jbv));
    for (i = 0; i < 16; i++) {
        result.data[i] = (p.data[i]+q.data[i])/2;
        result.data[i+16] = ((p.data[i^bit]-q.data[i^bit])/(i & bit ? even : odd))/2;
    }
    return result;
}

static inline multicomplex5
multicomplex5_pow(vec5 k, multicomplex5 a, multicomplex5 b)
{
    return multicomplex5_rpow(k, a, b);
}

static inline multicomplex5
multicomplex5_lpow(vec5 k, multicomplex5 a, multicomplex5 b)
{
    return multicomplex5_rpow(k, a, b);
}

static inline multicomplex5
multicomplex5_square(vec5 k, multicomplex5 a)
{
    return multicomplex5_mul(k, a, a);
}

static inline multicomplex5
multicomplex5_radd(double a, multicomplex5 b)
{
    return multicomplex5_add(multicomplex5_from_re(a), b);
}

static inline multicomplex5
multicomplex5_ladd(multicomplex5 a, double b)
{
    return multicomplex5_add(a, multicomplex5_from_re(b));
}

static inline multicomplex5
multicomplex5_rmul(vec5 k, double a, multicomplex5 b)
{
    return multicomplex5_mul(k, multicomplex5_from_re(a), b);
}

static inline multicomplex5
multicomplex5_lmul(vec5 k, multicomplex5 a, double b)
{
    return multicomplex5_mul(k, a, multicomplex5_from_re(b));
}

static inline multicomplex5
multicomplex5_rdiv(vec5 k, double a, multicomplex5 b)
{
    return multicomplex5_div(k, multicomplex5_from_re(a), b);
}

static inline multicomplex5
multicomplex5_ldiv(vec5 k, multicomplex5 a, double b)
{
    return multicomplex5_div(k, a, multicomplex5_from_re(b));
}

static inline multicomplex5
multicomplex5_rpows(vec5 k, double a, multicomplex5 b)
{
    return multicomplex5_rpow(k, multicomplex5_from_re(a), b);
}

static inline multicomplex5
multicomplex5_lpows(vec5 k, multicomplex5 a, double b)
{
    return multicomplex5_rpow(k, a, multicomplex5_from_re(b));
}

static inline multicomplex5
multicomplex5_expm1(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_expm1(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_expm1(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_expm1(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_expm1, multicomplex4_expm1_deriv1);
}

static inline multicomplex5
multicomplex5_log1p(vec5 k, multicomplex5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return complex5_log1p(a);
    }
    if (k.data[0] == 1 && k.data[1] == 1 && k.data[2] == 1 && k.data[3] == 1 && k.data[4] == 1) {
        return splitcomplex5_log1p(a);
    }
    if (k.data[0] == 0 && k.data[1] == 0 && k.data[2] == 0 && k.data[3] == 0 && k.data[4] == 0) {
        return dual5_log1p(a);
    }
    return multicomplex5_apply(k, a, multicomplex4_log1p, multicomplex4_log1p_deriv1);
}

static inline multicomplex5
multicomplex5_mix(vec5 k, multicomplex5 a, multicomplex5 b, double t)
{
    multicomplex5 delta = multicomplex5_mul(k, multicomplex5_inv(k, a), b);
    return multicomplex5_mul(k, a, multicomplex5_lpows(k, delta, t));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
multicomplex5_lmat(vec5 k, multicomplex5 a, double result[ALGPY_MULTICOMPLEX5_N*ALGPY_MULTICOMPLEX5_N])
{
    multicomplex5 basis = {{0}};
    multicomplex5 product;
    int row, col;
    for (col = 0; col < ALGPY_MULTICOMPLEX5_N; col++) {
        basis.data[col] = 1;
        product = multicomplex5_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_MULTICOMPLEX5_N; row++) {
            result[row * ALGPY_MULTICOMPLEX5_N + col] = product.data[row];
        }
    }
}

static inline void
multicomplex5_rmat(vec5 k, multicomplex5 a, double result[ALGPY_MULTICOMPLEX5_N*ALGPY_MULTICOMPLEX5_N])
{
    multicomplex5 basis = {{0}};
    multicomplex5 product;
    int row, col;
    for (row = 0; row < ALGPY_MULTICOMPLEX5_N; row++) {
        basis.data[row] = 1;
        product = multicomplex5_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_MULTICOMPLEX5_N; col++) {
            result[row * ALGPY_MULTICOMPLEX5_N + col] = product.data[col];
        }
    }
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline multicomplex5
multicomplex5_grade0(multicomplex5 a)
{
    const double *x = a.data;
    return multicomplex5_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline multicomplex5
multicomplex5_grade1(multicomplex5 a)
{
    const double *x = a.data;
    return multicomplex5_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline multicomplex5
multicomplex5_grade2(multicomplex5 a)
{
    const double *x = a.data;
    return multicomplex5_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline multicomplex5
multicomplex5_grade3(multicomplex5 a)
{
    const double *x = a.data;
    return multicomplex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline multicomplex5
multicomplex5_grade4(multicomplex5 a)
{
    const double *x = a.data;
    return multicomplex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline multicomplex5
multicomplex5_grade5(multicomplex5 a)
{
    const double *x = a.data;
    return multicomplex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

#endif /* ALGPY_MULTICOMPLEX5_H */
