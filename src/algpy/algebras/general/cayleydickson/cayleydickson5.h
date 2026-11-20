/*
 * Core arithmetic of the Cayley-Dickson algebra with five generators.
 * The parameters k are passed explicitly to every operation that depends on
 * them. This header has no Python or NumPy dependencies.
 */
#ifndef ALGPY_CAYLEYDICKSON5_H
#define ALGPY_CAYLEYDICKSON5_H

#define ALGPY_CAYLEYDICKSON5_N 32

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../../fields/complex.h"
#include "../../../fields/double.h"
#include "../../cayleydickson/trigintaduonion.h"
#include "cayleydickson4.h"
#include "../../../vectorspaces/vec5.h"
#include "../../../vectorspaces/vec32.h"

typedef vec32 cayleydickson5;

static inline cayleydickson5
cayleydickson5_make(const double coeffs[ALGPY_CAYLEYDICKSON5_N])
{
    cayleydickson5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline cayleydickson5
cayleydickson5_re(cayleydickson5 a)
{
    cayleydickson5 result = {{0}};
    result.data[0] = a.data[0];
    return result;
}

static inline cayleydickson5
cayleydickson5_im(cayleydickson5 a)
{
    return cayleydickson5_make((const double[]){0, a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15], a.data[16], a.data[17], a.data[18], a.data[19], a.data[20], a.data[21], a.data[22], a.data[23], a.data[24], a.data[25], a.data[26], a.data[27], a.data[28], a.data[29], a.data[30], a.data[31]});
}

static inline cayleydickson5
cayleydickson5_from_re(double a)
{
    cayleydickson5 result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
cayleydickson5_to_re(cayleydickson5 a)
{
    return a.data[0];
}

static inline cayleydickson5
cayleydickson5_one(void)
{
    return cayleydickson5_from_re(1);
}

/* Cayley-Dickson doubling of the lower four-generator algebra. */
static inline cayleydickson5
cayleydickson5_mul(vec5 k, cayleydickson5 a, cayleydickson5 b)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return trigintaduonion_mul(a, b);
    }
    vec4 base = {{k.data[0], k.data[1], k.data[2], k.data[3]}};
    cayleydickson4 x0 = {{a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]}};
    cayleydickson4 x1 = {{a.data[16], a.data[17], a.data[18], a.data[19], a.data[20], a.data[21], a.data[22], a.data[23], a.data[24], a.data[25], a.data[26], a.data[27], a.data[28], a.data[29], a.data[30], a.data[31]}};
    cayleydickson4 y0 = {{b.data[0], b.data[1], b.data[2], b.data[3], b.data[4], b.data[5], b.data[6], b.data[7], b.data[8], b.data[9], b.data[10], b.data[11], b.data[12], b.data[13], b.data[14], b.data[15]}};
    cayleydickson4 y1 = {{b.data[16], b.data[17], b.data[18], b.data[19], b.data[20], b.data[21], b.data[22], b.data[23], b.data[24], b.data[25], b.data[26], b.data[27], b.data[28], b.data[29], b.data[30], b.data[31]}};
    cayleydickson4 left = vec16_add(
        cayleydickson4_mul(base, x0, y0),
        vec16_scale(cayleydickson4_mul(base, cayleydickson4_conjugate(y1), x1), k.data[4]));
    cayleydickson4 right = vec16_add(
        cayleydickson4_mul(base, y1, x0),
        cayleydickson4_mul(base, x1, cayleydickson4_conjugate(y0)));
    cayleydickson5 result;
    memcpy(result.data, left.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_CAYLEYDICKSON5_N/2, right.data, sizeof(result.data)/2);
    return result;
}

static inline cayleydickson5
cayleydickson5_conjugate(cayleydickson5 a)
{
    return vec32_sub(cayleydickson5_re(a), cayleydickson5_im(a));
}

/* norm2(a) = sum_S product_{i in S}(-k[i]) * a[S]^2. */
static inline double
cayleydickson5_norm2(vec5 k, cayleydickson5 a)
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
        + k.data[0]*k.data[1]*k.data[2]*k.data[3]*a.data[15]*a.data[15]
        - k.data[4]*a.data[16]*a.data[16]
        + k.data[0]*k.data[4]*a.data[17]*a.data[17]
        + k.data[1]*k.data[4]*a.data[18]*a.data[18]
        - k.data[0]*k.data[1]*k.data[4]*a.data[19]*a.data[19]
        + k.data[2]*k.data[4]*a.data[20]*a.data[20]
        - k.data[0]*k.data[2]*k.data[4]*a.data[21]*a.data[21]
        - k.data[1]*k.data[2]*k.data[4]*a.data[22]*a.data[22]
        + k.data[0]*k.data[1]*k.data[2]*k.data[4]*a.data[23]*a.data[23]
        + k.data[3]*k.data[4]*a.data[24]*a.data[24]
        - k.data[0]*k.data[3]*k.data[4]*a.data[25]*a.data[25]
        - k.data[1]*k.data[3]*k.data[4]*a.data[26]*a.data[26]
        + k.data[0]*k.data[1]*k.data[3]*k.data[4]*a.data[27]*a.data[27]
        - k.data[2]*k.data[3]*k.data[4]*a.data[28]*a.data[28]
        + k.data[0]*k.data[2]*k.data[3]*k.data[4]*a.data[29]*a.data[29]
        + k.data[1]*k.data[2]*k.data[3]*k.data[4]*a.data[30]*a.data[30]
        - k.data[0]*k.data[1]*k.data[2]*k.data[3]*k.data[4]*a.data[31]*a.data[31];
}

static inline cayleydickson5
cayleydickson5_inv(vec5 k, cayleydickson5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return trigintaduonion_inv(a);
    }
    return vec32_div_scalar(cayleydickson5_conjugate(a),
                            cayleydickson5_norm2(k, a));
}

static inline cayleydickson5
cayleydickson5_div(vec5 k, cayleydickson5 a, cayleydickson5 b)
{
    return cayleydickson5_mul(k, a, cayleydickson5_inv(k, b));
}

static inline void
cayleydickson5_basis(cayleydickson5 result[ALGPY_CAYLEYDICKSON5_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON5_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON5_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson5_gens(cayleydickson5 result[5])
{
    result[0] = cayleydickson5_make((const double[ALGPY_CAYLEYDICKSON5_N]){[1] = 1});
    result[1] = cayleydickson5_make((const double[ALGPY_CAYLEYDICKSON5_N]){[2] = 1});
    result[2] = cayleydickson5_make((const double[ALGPY_CAYLEYDICKSON5_N]){[4] = 1});
    result[3] = cayleydickson5_make((const double[ALGPY_CAYLEYDICKSON5_N]){[8] = 1});
    result[4] = cayleydickson5_make((const double[ALGPY_CAYLEYDICKSON5_N]){[16] = 1});
}

static inline void
cayleydickson5_coeffs(
    vec5 k,
    cayleydickson5 result[ALGPY_CAYLEYDICKSON5_N *
                          ALGPY_CAYLEYDICKSON5_N])
{
    cayleydickson5 basis[ALGPY_CAYLEYDICKSON5_N];
    int i, j;
    cayleydickson5_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON5_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON5_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON5_N + j] =
                cayleydickson5_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
cayleydickson5_minpoly(vec5 k, cayleydickson5 a, double result[3])
{
    result[0] = cayleydickson5_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

static inline void
cayleydickson5_eivals(vec5 k, cayleydickson5 a, complex result[2])
{
    complex s = complex_sqrt_real(-cayleydickson5_norm2(k, cayleydickson5_im(a)));
    result[0] = complex_make((const double[]){a.data[0] - s.data[0], -s.data[1]});
    result[1] = complex_make((const double[]){a.data[0] + s.data[0], +s.data[1]});
}

/* Common algebra API. */
static inline cayleydickson5
cayleydickson5_zero(void)
{
    return vec32_zero();
}

static inline int
cayleydickson5_is_real(cayleydickson5 a)
{
    return vec32_is_real(a);
}

static inline cayleydickson5
cayleydickson5_pos(cayleydickson5 a)
{
    return vec32_pos(a);
}

static inline cayleydickson5
cayleydickson5_neg(cayleydickson5 a)
{
    return vec32_neg(a);
}

static inline cayleydickson5
cayleydickson5_add(cayleydickson5 a, cayleydickson5 b)
{
    return vec32_add(a, b);
}

static inline cayleydickson5
cayleydickson5_sub(cayleydickson5 a, cayleydickson5 b)
{
    return vec32_sub(a, b);
}

static inline int
cayleydickson5_equal(cayleydickson5 a, cayleydickson5 b)
{
    return vec32_all_eq(a, b);
}

static inline int
cayleydickson5_not_equal(cayleydickson5 a, cayleydickson5 b)
{
    return !cayleydickson5_equal(a, b);
}

/* Every element satisfies a^2 - 2*re(a)*a + norm2(a) = 0.
 * For a = w + v this gives v^2 = q: complex if q < 0, split if q > 0,
 * and f(a) = f(w) + f'(w)*v if q == 0. */
static inline cayleydickson5
cayleydickson5_slice(cayleydickson5 a, double real, double scale)
{
    cayleydickson5 result;
    int i;
    result.data[0] = real;
    for (i = 1; i < ALGPY_CAYLEYDICKSON5_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : scale*a.data[i];
    }
    return result;
}

static inline cayleydickson5
cayleydickson5_real_slice(vec5 k, complex value)
{
    cayleydickson5 result = cayleydickson5_from_re(value.data[0]);
    int i, p;
    if (value.data[1] == 0) {
        return result;
    }
    for (i = 1; i < ALGPY_CAYLEYDICKSON5_N; i++) {
        double square = -1;
        for (p = 0; p < 5; p++) {
            if (i & (1 << p)) {
                square *= -k.data[p];
            }
        }
        if (square < 0) {
            result.data[i] = value.data[1]/sqrt(-square);
            return result;
        }
    }
    return cayleydickson5_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline cayleydickson5
cayleydickson5_apply(vec5 k, cayleydickson5 a,
    double (*f)(double), double (*deriv1)(double), complex (*cf)(complex))
{
    double w = a.data[0], q, radius, plus, minus;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return trigintaduonion_apply(a, cf);
    }
    if (vec32_is_real(a)) {
        return cayleydickson5_real_slice(k, cf(complex_make((const double[]){w, 0})));
    }
    q = -cayleydickson5_norm2(k, cayleydickson5_im(a));
    if (q < 0) {
        complex value;
        cayleydickson5 result;
        radius = sqrt(-q);
        value = cf(complex_make((const double[]){w, radius}));
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON5_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = f(w + radius);
        minus = f(w - radius);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson5_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
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
                return cayleydickson5_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson5_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        double real = f(w), scale = deriv1(w);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson5_slice(a, real, scale);
        }
    }
    return cayleydickson5_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline cayleydickson5
cayleydickson5_exp(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_exp, double_exp_deriv1, complex_exp);
}

static inline cayleydickson5
cayleydickson5_log(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_log, double_log_deriv1, complex_log);
}

static inline cayleydickson5
cayleydickson5_sqrt(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_sqrt, double_sqrt_deriv1, complex_sqrt);
}

static inline cayleydickson5
cayleydickson5_cos(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_cos, double_cos_deriv1, complex_cos);
}

static inline cayleydickson5
cayleydickson5_sin(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_sin, double_sin_deriv1, complex_sin);
}

static inline cayleydickson5
cayleydickson5_tan(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_tan, double_tan_deriv1, complex_tan);
}

static inline cayleydickson5
cayleydickson5_acos(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_acos, double_acos_deriv1, complex_acos);
}

static inline cayleydickson5
cayleydickson5_asin(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_asin, double_asin_deriv1, complex_asin);
}

static inline cayleydickson5
cayleydickson5_atan(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_atan, double_atan_deriv1, complex_atan);
}

static inline cayleydickson5
cayleydickson5_cosh(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_cosh, double_cosh_deriv1, complex_cosh);
}

static inline cayleydickson5
cayleydickson5_sinh(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_sinh, double_sinh_deriv1, complex_sinh);
}

static inline cayleydickson5
cayleydickson5_tanh(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_tanh, double_tanh_deriv1, complex_tanh);
}

static inline cayleydickson5
cayleydickson5_acosh(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_acosh, double_acosh_deriv1, complex_acosh);
}

static inline cayleydickson5
cayleydickson5_asinh(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_asinh, double_asinh_deriv1, complex_asinh);
}

static inline cayleydickson5
cayleydickson5_atanh(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_atanh, double_atanh_deriv1, complex_atanh);
}

static inline cayleydickson5
cayleydickson5_square(vec5 k, cayleydickson5 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return trigintaduonion_square(a);
    }
    return cayleydickson5_slice(a, a.data[0]*a.data[0] - cayleydickson5_norm2(k, cayleydickson5_im(a)), 2*a.data[0]);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline cayleydickson5
cayleydickson5_ipow(vec5 k, cayleydickson5 a, int64_t b)
{
    cayleydickson5 result = cayleydickson5_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson5_mul(k, a, a);
    }
    if (b == -1) {
        return cayleydickson5_inv(k, a);
    }
    if (b < 0) {
        a = cayleydickson5_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = cayleydickson5_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = cayleydickson5_mul(k, a, a);
        }
    }
    return result;
}

static inline cayleydickson5
cayleydickson5_lpows(vec5 k, cayleydickson5 a, double b)
{
    double w = a.data[0], q, radius, plus, minus, real, scale;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1 && k.data[3] == -1 && k.data[4] == -1) {
        return trigintaduonion_lpows(a, b);
    }
    if (b == 0) {
        return cayleydickson5_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson5_square(k, a);
    }
    if (vec32_is_real(a)) {
        real = pow(w, b);
        if (!isnan(real)) {
            return cayleydickson5_from_re(real);
        }
        return cayleydickson5_real_slice(k, complex_lpows(complex_make((const double[]){w, 0}), b));
    }
    q = -cayleydickson5_norm2(k, cayleydickson5_im(a));
    if (q < 0) {
        complex value;
        cayleydickson5 result;
        radius = sqrt(-q);
        value = complex_lpows(complex_make((const double[]){w, radius}), b);
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON5_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = pow(w + radius, b);
        minus = pow(w - radius, b);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson5_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
        }
        if (radius < 1e-8*fabs(w)/fmax(1, fabs(b))) {
            real = pow(w, b);
            scale = b*pow(w, b - 1);
            if (isfinite(real) && isfinite(scale)) {
                return cayleydickson5_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson5_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        real = pow(w, b);
        scale = b*pow(w, b - 1);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson5_slice(a, real, scale);
        }
    }
    return cayleydickson5_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline cayleydickson5
cayleydickson5_rpow(vec5 k, cayleydickson5 a, cayleydickson5 b)
{
    if (vec32_is_real(b)) {
        return cayleydickson5_lpows(k, a, b.data[0]);
    }
    return cayleydickson5_exp(k, cayleydickson5_mul(k, cayleydickson5_log(k, a), b));
}

static inline cayleydickson5
cayleydickson5_pow(vec5 k, cayleydickson5 a, cayleydickson5 b)
{
    return cayleydickson5_rpow(k, a, b);
}

static inline cayleydickson5
cayleydickson5_lpow(vec5 k, cayleydickson5 a, cayleydickson5 b)
{
    return cayleydickson5_exp(k, cayleydickson5_mul(k, b, cayleydickson5_log(k, a)));
}

static inline cayleydickson5
cayleydickson5_radd(double a, cayleydickson5 b)
{
    return cayleydickson5_add(cayleydickson5_from_re(a), b);
}

static inline cayleydickson5
cayleydickson5_ladd(cayleydickson5 a, double b)
{
    return cayleydickson5_add(a, cayleydickson5_from_re(b));
}

static inline cayleydickson5
cayleydickson5_rsub(double a, cayleydickson5 b)
{
    return cayleydickson5_sub(cayleydickson5_from_re(a), b);
}

static inline cayleydickson5
cayleydickson5_lsub(cayleydickson5 a, double b)
{
    return cayleydickson5_sub(a, cayleydickson5_from_re(b));
}

static inline cayleydickson5
cayleydickson5_rmul(vec5 k, double a, cayleydickson5 b)
{
    return cayleydickson5_mul(k, cayleydickson5_from_re(a), b);
}

static inline cayleydickson5
cayleydickson5_lmul(vec5 k, cayleydickson5 a, double b)
{
    return cayleydickson5_mul(k, a, cayleydickson5_from_re(b));
}

static inline cayleydickson5
cayleydickson5_rdiv(vec5 k, double a, cayleydickson5 b)
{
    return cayleydickson5_div(k, cayleydickson5_from_re(a), b);
}

static inline cayleydickson5
cayleydickson5_ldiv(vec5 k, cayleydickson5 a, double b)
{
    return cayleydickson5_div(k, a, cayleydickson5_from_re(b));
}

static inline cayleydickson5
cayleydickson5_rpows(vec5 k, double a, cayleydickson5 b)
{
    return cayleydickson5_rpow(k, cayleydickson5_from_re(a), b);
}

static inline cayleydickson5
cayleydickson5_expm1(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_expm1, double_expm1_deriv1, complex_expm1);
}

static inline cayleydickson5
cayleydickson5_log1p(vec5 k, cayleydickson5 a)
{
    return cayleydickson5_apply(k, a, double_log1p, double_log1p_deriv1, complex_log1p);
}

static inline cayleydickson5
cayleydickson5_mix(vec5 k, cayleydickson5 a, cayleydickson5 b, double t)
{
    cayleydickson5 delta = cayleydickson5_mul(k, cayleydickson5_inv(k, a), b);
    return cayleydickson5_mul(k, a, cayleydickson5_lpows(k, delta, t));
}

static inline complex
cayleydickson5_norm(vec5 k, cayleydickson5 a)
{
    return complex_root_real(cayleydickson5_norm2(k, a), 2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
cayleydickson5_lmat(vec5 k, cayleydickson5 a, double result[ALGPY_CAYLEYDICKSON5_N*ALGPY_CAYLEYDICKSON5_N])
{
    cayleydickson5 basis = {{0}};
    cayleydickson5 product;
    int row, col;
    for (col = 0; col < ALGPY_CAYLEYDICKSON5_N; col++) {
        basis.data[col] = 1;
        product = cayleydickson5_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CAYLEYDICKSON5_N; row++) {
            result[row * ALGPY_CAYLEYDICKSON5_N + col] = product.data[row];
        }
    }
}

static inline void
cayleydickson5_rmat(vec5 k, cayleydickson5 a, double result[ALGPY_CAYLEYDICKSON5_N*ALGPY_CAYLEYDICKSON5_N])
{
    cayleydickson5 basis = {{0}};
    cayleydickson5 product;
    int row, col;
    for (row = 0; row < ALGPY_CAYLEYDICKSON5_N; row++) {
        basis.data[row] = 1;
        product = cayleydickson5_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CAYLEYDICKSON5_N; col++) {
            result[row * ALGPY_CAYLEYDICKSON5_N + col] = product.data[col];
        }
    }
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline cayleydickson5
cayleydickson5_grade0(cayleydickson5 a)
{
    const double *x = a.data;
    return cayleydickson5_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline cayleydickson5
cayleydickson5_grade1(cayleydickson5 a)
{
    const double *x = a.data;
    return cayleydickson5_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline cayleydickson5
cayleydickson5_grade2(cayleydickson5 a)
{
    const double *x = a.data;
    return cayleydickson5_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline cayleydickson5
cayleydickson5_grade3(cayleydickson5 a)
{
    const double *x = a.data;
    return cayleydickson5_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline cayleydickson5
cayleydickson5_grade4(cayleydickson5 a)
{
    const double *x = a.data;
    return cayleydickson5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline cayleydickson5
cayleydickson5_grade5(cayleydickson5 a)
{
    const double *x = a.data;
    return cayleydickson5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

#endif /* ALGPY_CAYLEYDICKSON5_H */
