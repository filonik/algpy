/*
 * Core arithmetic of the Cayley-Dickson algebra with four generators.
 * The parameters k are passed explicitly to every operation that depends on
 * them. This header has no Python or NumPy dependencies.
 */
#ifndef ALGPY_CAYLEYDICKSON4_H
#define ALGPY_CAYLEYDICKSON4_H

#define ALGPY_CAYLEYDICKSON4_N 16

#include <stdint.h>
#include <math.h>

#include "../../../fields/complex.h"
#include "../../../fields/double.h"
#include "../../cayleydickson/sedenion.h"
#include "cayleydickson3.h"
#include "../../../vectorspaces/vec4.h"
#include "../../../vectorspaces/vec16.h"

typedef vec16 cayleydickson4;

static inline cayleydickson4
cayleydickson4_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012,
                    double x3, double x03, double x13, double x013,
                    double x23, double x023, double x123, double x0123)
{
    cayleydickson4 r = {{x, x0, x1, x01, x2, x02, x12, x012,
                         x3, x03, x13, x013, x23, x023, x123, x0123}};
    return r;
}

static inline cayleydickson4
cayleydickson4_re(cayleydickson4 a)
{
    return cayleydickson4_make(a.data[0], 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0);
}

static inline cayleydickson4
cayleydickson4_im(cayleydickson4 a)
{
    return cayleydickson4_make(0, a.data[1], a.data[2], a.data[3],
                               a.data[4], a.data[5], a.data[6], a.data[7],
                               a.data[8], a.data[9], a.data[10], a.data[11],
                               a.data[12], a.data[13], a.data[14], a.data[15]);
}

static inline cayleydickson4
cayleydickson4_from_re(double a)
{
    return cayleydickson4_make(a, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
cayleydickson4_to_re(cayleydickson4 a)
{
    return a.data[0];
}

static inline cayleydickson4
cayleydickson4_one(void)
{
    return cayleydickson4_from_re(1);
}

static inline cayleydickson4
cayleydickson4_mul(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return sedenion_mul(a, b);
    }
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    cayleydickson3 x0 = {{a.data[0], a.data[1], a.data[2], a.data[3],
                          a.data[4], a.data[5], a.data[6], a.data[7]}};
    cayleydickson3 x1 = {{a.data[8], a.data[9], a.data[10], a.data[11],
                          a.data[12], a.data[13], a.data[14], a.data[15]}};
    cayleydickson3 y0 = {{b.data[0], b.data[1], b.data[2], b.data[3],
                          b.data[4], b.data[5], b.data[6], b.data[7]}};
    cayleydickson3 y1 = {{b.data[8], b.data[9], b.data[10], b.data[11],
                          b.data[12], b.data[13], b.data[14], b.data[15]}};
    cayleydickson3 left = vec8_add(
        cayleydickson3_mul(k012, x0, y0),
        vec8_scale(cayleydickson3_mul(
                       k012, cayleydickson3_conjugate(y1), x1),
                   k.data[3]));
    cayleydickson3 right = vec8_add(
        cayleydickson3_mul(k012, y1, x0),
        cayleydickson3_mul(k012, x1, cayleydickson3_conjugate(y0)));
    return cayleydickson4_make(
        left.data[0], left.data[1], left.data[2], left.data[3],
        left.data[4], left.data[5], left.data[6], left.data[7],
        right.data[0], right.data[1], right.data[2], right.data[3],
        right.data[4], right.data[5], right.data[6], right.data[7]);
}

static inline cayleydickson4
cayleydickson4_conjugate(cayleydickson4 a)
{
    return vec16_sub(cayleydickson4_re(a), cayleydickson4_im(a));
}

static inline double
cayleydickson4_norm2(vec4 k, cayleydickson4 a)
{
    return a.data[0]*a.data[0]
        - k.data[0]*a.data[1]*a.data[1]
        - k.data[1]*a.data[2]*a.data[2]
        + k.data[0]*k.data[1]*a.data[3]*a.data[3]
        - k.data[2]*a.data[4]*a.data[4]
        + k.data[0]*k.data[2]*a.data[5]*a.data[5]
        + k.data[1]*k.data[2]*a.data[6]*a.data[6]
        - k.data[0]*k.data[1]*k.data[2]*a.data[7]*a.data[7]
        - k.data[3]*a.data[8]*a.data[8]
        + k.data[0]*k.data[3]*a.data[9]*a.data[9]
        + k.data[1]*k.data[3]*a.data[10]*a.data[10]
        - k.data[0]*k.data[1]*k.data[3]*a.data[11]*a.data[11]
        + k.data[2]*k.data[3]*a.data[12]*a.data[12]
        - k.data[0]*k.data[2]*k.data[3]*a.data[13]*a.data[13]
        - k.data[1]*k.data[2]*k.data[3]*a.data[14]*a.data[14]
        + k.data[0]*k.data[1]*k.data[2]*k.data[3]*a.data[15]*a.data[15];
}

static inline cayleydickson4
cayleydickson4_inv(vec4 k, cayleydickson4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return sedenion_inv(a);
    }
    return vec16_div_scalar(cayleydickson4_conjugate(a),
                            cayleydickson4_norm2(k, a));
}

static inline cayleydickson4
cayleydickson4_div(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    return cayleydickson4_mul(k, a, cayleydickson4_inv(k, b));
}

static inline void
cayleydickson4_basis(cayleydickson4 result[ALGPY_CAYLEYDICKSON4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson4_coeffs(
    vec4 k,
    cayleydickson4 result[ALGPY_CAYLEYDICKSON4_N *
                          ALGPY_CAYLEYDICKSON4_N])
{
    cayleydickson4 basis[ALGPY_CAYLEYDICKSON4_N];
    int i, j;
    cayleydickson4_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON4_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON4_N + j] =
                cayleydickson4_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
cayleydickson4_minpoly(vec4 k, cayleydickson4 a, double result[3])
{
    result[0] = cayleydickson4_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

static inline void
cayleydickson4_eivals(vec4 k, cayleydickson4 a, complex result[2])
{
    complex s = complex_sqrt_real(-cayleydickson4_norm2(k, cayleydickson4_im(a)));
    result[0] = complex_make(a.data[0] - s.data[0], -s.data[1]);
    result[1] = complex_make(a.data[0] + s.data[0], +s.data[1]);
}

/* Common algebra API. */
static inline cayleydickson4
cayleydickson4_zero(void)
{
    return vec16_zero();
}

static inline int
cayleydickson4_is_real(cayleydickson4 a)
{
    return vec16_is_real(a);
}

static inline cayleydickson4
cayleydickson4_pos(cayleydickson4 a)
{
    return vec16_pos(a);
}

static inline cayleydickson4
cayleydickson4_neg(cayleydickson4 a)
{
    return vec16_neg(a);
}

static inline cayleydickson4
cayleydickson4_add(cayleydickson4 a, cayleydickson4 b)
{
    return vec16_add(a, b);
}

static inline cayleydickson4
cayleydickson4_sub(cayleydickson4 a, cayleydickson4 b)
{
    return vec16_sub(a, b);
}

static inline int
cayleydickson4_equal(cayleydickson4 a, cayleydickson4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
cayleydickson4_not_equal(cayleydickson4 a, cayleydickson4 b)
{
    return !cayleydickson4_equal(a, b);
}

/* Every element satisfies a^2 - 2*re(a)*a + norm2(a) = 0.
 * For a = w + v this gives v^2 = q: complex if q < 0, split if q > 0,
 * and f(a) = f(w) + f'(w)*v if q == 0. */
static inline cayleydickson4
cayleydickson4_slice(cayleydickson4 a, double real, double scale)
{
    cayleydickson4 result;
    int i;
    result.data[0] = real;
    for (i = 1; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : scale*a.data[i];
    }
    return result;
}

static inline cayleydickson4
cayleydickson4_real_slice(vec4 k, complex value)
{
    cayleydickson4 result = cayleydickson4_from_re(value.data[0]);
    int i, p;
    if (value.data[1] == 0) {
        return result;
    }
    for (i = 1; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        double square = -1;
        for (p = 0; p < 4; p++) {
            if (i & (1 << p)) {
                square *= -k.data[p];
            }
        }
        if (square < 0) {
            result.data[i] = value.data[1]/sqrt(-square);
            return result;
        }
    }
    return cayleydickson4_make(
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_apply(vec4 k, cayleydickson4 a,
    double (*f)(double), double (*deriv1)(double), complex (*cf)(complex))
{
    double w = a.data[0], q, radius, plus, minus;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return sedenion_apply(a, cf);
    }
    if (vec16_is_real(a)) {
        return cayleydickson4_real_slice(k, cf(complex_make(w, 0)));
    }
    q = -cayleydickson4_norm2(k, cayleydickson4_im(a));
    if (q < 0) {
        complex value;
        cayleydickson4 result;
        radius = sqrt(-q);
        value = cf(complex_make(w, radius));
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON4_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = f(w + radius);
        minus = f(w - radius);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
        }
        /* Native real functions have a smooth derivative limit when the
         * roots nearly coincide. Respect their nearest branch point or pole. */
        double distance = 1;
        if (f == double_log || f == double_sqrt) {
            distance = fabs(w);
        } else if (f == double_log1p) {
            distance = fabs(1 + w);
        } else if (f == double_asin || f == double_acos ||
                   f == double_acosh || f == double_atanh) {
            distance = fmin(fabs(1 - w), fabs(1 + w));
        } else if (f == double_atan || f == double_asinh) {
            distance = hypot(1, w);
        } else if (f == double_tan) {
            distance = fabs(remainder(w - M_PI/2, M_PI));
        }
        if (radius < 1e-8*distance &&
            (f == double_exp ||
             f == double_log ||
             f == double_sqrt ||
             f == double_cos ||
             f == double_sin ||
             f == double_tan ||
             f == double_acos ||
             f == double_asin ||
             f == double_atan ||
             f == double_cosh ||
             f == double_sinh ||
             f == double_tanh ||
             f == double_acosh ||
             f == double_asinh ||
             f == double_atanh ||
             f == double_expm1 ||
             f == double_log1p)) {
            double real = f(w), scale = deriv1(w);
            if (isfinite(real) && isfinite(scale)) {
                return cayleydickson4_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson4_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        double real = f(w), scale = deriv1(w);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson4_slice(a, real, scale);
        }
    }
    return cayleydickson4_make(
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_exp(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_exp, double_exp_deriv1, complex_exp);
}

static inline cayleydickson4
cayleydickson4_log(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_log, double_log_deriv1, complex_log);
}

static inline cayleydickson4
cayleydickson4_sqrt(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_sqrt, double_sqrt_deriv1, complex_sqrt);
}

static inline cayleydickson4
cayleydickson4_cos(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_cos, double_cos_deriv1, complex_cos);
}

static inline cayleydickson4
cayleydickson4_sin(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_sin, double_sin_deriv1, complex_sin);
}

static inline cayleydickson4
cayleydickson4_tan(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_tan, double_tan_deriv1, complex_tan);
}

static inline cayleydickson4
cayleydickson4_acos(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_acos, double_acos_deriv1, complex_acos);
}

static inline cayleydickson4
cayleydickson4_asin(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_asin, double_asin_deriv1, complex_asin);
}

static inline cayleydickson4
cayleydickson4_atan(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_atan, double_atan_deriv1, complex_atan);
}

static inline cayleydickson4
cayleydickson4_cosh(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_cosh, double_cosh_deriv1, complex_cosh);
}

static inline cayleydickson4
cayleydickson4_sinh(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_sinh, double_sinh_deriv1, complex_sinh);
}

static inline cayleydickson4
cayleydickson4_tanh(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_tanh, double_tanh_deriv1, complex_tanh);
}

static inline cayleydickson4
cayleydickson4_acosh(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_acosh, double_acosh_deriv1, complex_acosh);
}

static inline cayleydickson4
cayleydickson4_asinh(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_asinh, double_asinh_deriv1, complex_asinh);
}

static inline cayleydickson4
cayleydickson4_atanh(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_atanh, double_atanh_deriv1, complex_atanh);
}

static inline cayleydickson4
cayleydickson4_square(vec4 k, cayleydickson4 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return sedenion_square(a);
    }
    return cayleydickson4_slice(a, a.data[0]*a.data[0] - cayleydickson4_norm2(k, cayleydickson4_im(a)), 2*a.data[0]);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline cayleydickson4
cayleydickson4_ipow(vec4 k, cayleydickson4 a, int64_t b)
{
    cayleydickson4 result = cayleydickson4_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson4_mul(k, a, a);
    }
    if (b == -1) {
        return cayleydickson4_inv(k, a);
    }
    if (b < 0) {
        a = cayleydickson4_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = cayleydickson4_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = cayleydickson4_mul(k, a, a);
        }
    }
    return result;
}

static inline cayleydickson4
cayleydickson4_lpows(vec4 k, cayleydickson4 a, double b)
{
    double w = a.data[0], q, radius, plus, minus, real, scale;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1) {
        return sedenion_lpows(a, b);
    }
    if (b == 0) {
        return cayleydickson4_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson4_square(k, a);
    }
    if (vec16_is_real(a)) {
        real = pow(w, b);
        if (!isnan(real)) {
            return cayleydickson4_from_re(real);
        }
        return cayleydickson4_real_slice(k, complex_lpows(complex_make(w, 0), b));
    }
    q = -cayleydickson4_norm2(k, cayleydickson4_im(a));
    if (q < 0) {
        complex value;
        cayleydickson4 result;
        radius = sqrt(-q);
        value = complex_lpows(complex_make(w, radius), b);
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON4_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = pow(w + radius, b);
        minus = pow(w - radius, b);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
        }
        if (radius < 1e-8*fabs(w)/fmax(1, fabs(b))) {
            real = pow(w, b);
            scale = b*pow(w, b - 1);
            if (isfinite(real) && isfinite(scale)) {
                return cayleydickson4_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson4_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        real = pow(w, b);
        scale = b*pow(w, b - 1);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson4_slice(a, real, scale);
        }
    }
    return cayleydickson4_make(
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN,
        NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_rpow(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    if (vec16_is_real(b)) {
        return cayleydickson4_lpows(k, a, b.data[0]);
    }
    return cayleydickson4_exp(k, cayleydickson4_mul(k, cayleydickson4_log(k, a), b));
}

static inline cayleydickson4
cayleydickson4_pow(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    return cayleydickson4_rpow(k, a, b);
}

static inline cayleydickson4
cayleydickson4_lpow(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    return cayleydickson4_exp(k, cayleydickson4_mul(k, b, cayleydickson4_log(k, a)));
}



static inline cayleydickson4
cayleydickson4_radd(double a, cayleydickson4 b)
{
    return cayleydickson4_add(cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_ladd(cayleydickson4 a, double b)
{
    return cayleydickson4_add(a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rsub(double a, cayleydickson4 b)
{
    return cayleydickson4_sub(cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_lsub(cayleydickson4 a, double b)
{
    return cayleydickson4_sub(a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rmul(vec4 k, double a, cayleydickson4 b)
{
    return cayleydickson4_mul(k, cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_lmul(vec4 k, cayleydickson4 a, double b)
{
    return cayleydickson4_mul(k, a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rdiv(vec4 k, double a, cayleydickson4 b)
{
    return cayleydickson4_div(k, cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_ldiv(vec4 k, cayleydickson4 a, double b)
{
    return cayleydickson4_div(k, a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rpows(vec4 k, double a, cayleydickson4 b)
{
    return cayleydickson4_rpow(k, cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_expm1(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_expm1, double_expm1_deriv1, complex_expm1);
}

static inline cayleydickson4
cayleydickson4_log1p(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_apply(k, a, double_log1p, double_log1p_deriv1, complex_log1p);
}

static inline cayleydickson4
cayleydickson4_mix(vec4 k, cayleydickson4 a, cayleydickson4 b, double t)
{
    cayleydickson4 delta = cayleydickson4_mul(k, cayleydickson4_inv(k, a), b);
    return cayleydickson4_mul(k, a, cayleydickson4_lpows(k, delta, t));
}

static inline complex
cayleydickson4_norm(vec4 k, cayleydickson4 a)
{
    return complex_root_real(cayleydickson4_norm2(k, a), 2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
cayleydickson4_lmat(vec4 k, cayleydickson4 a, double result[ALGPY_CAYLEYDICKSON4_N*ALGPY_CAYLEYDICKSON4_N])
{
    cayleydickson4 basis = {{0}};
    cayleydickson4 product;
    int row, col;
    for (col = 0; col < ALGPY_CAYLEYDICKSON4_N; col++) {
        basis.data[col] = 1;
        product = cayleydickson4_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CAYLEYDICKSON4_N; row++) {
            result[row * ALGPY_CAYLEYDICKSON4_N + col] = product.data[row];
        }
    }
}

static inline void
cayleydickson4_rmat(vec4 k, cayleydickson4 a, double result[ALGPY_CAYLEYDICKSON4_N*ALGPY_CAYLEYDICKSON4_N])
{
    cayleydickson4 basis = {{0}};
    cayleydickson4 product;
    int row, col;
    for (row = 0; row < ALGPY_CAYLEYDICKSON4_N; row++) {
        basis.data[row] = 1;
        product = cayleydickson4_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CAYLEYDICKSON4_N; col++) {
            result[row * ALGPY_CAYLEYDICKSON4_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_CAYLEYDICKSON4_H */
