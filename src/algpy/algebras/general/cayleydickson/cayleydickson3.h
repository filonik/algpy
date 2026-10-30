/*
 * Core arithmetic of the Cayley-Dickson algebra with three generators.
 * The parameters k are passed explicitly to every operation that depends on
 * them. This header has no Python or NumPy dependencies.
 */
#ifndef ALGPY_CAYLEYDICKSON3_H
#define ALGPY_CAYLEYDICKSON3_H

#define ALGPY_CAYLEYDICKSON3_N 8

#include <stdint.h>
#include <math.h>

#include "../../../fields/complex.h"
#include "../../../fields/double.h"
#include "../../cayleydickson/octonion.h"
#include "cayleydickson2.h"
#include "../../../vectorspaces/vec3.h"
#include "../../../vectorspaces/vec8.h"

typedef vec8 cayleydickson3;

static inline cayleydickson3
cayleydickson3_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012)
{
    cayleydickson3 r = {{x, x0, x1, x01, x2, x02, x12, x012}};
    return r;
}

static inline cayleydickson3
cayleydickson3_re(cayleydickson3 a)
{
    return cayleydickson3_make(a.data[0], 0, 0, 0, 0, 0, 0, 0);
}

static inline cayleydickson3
cayleydickson3_im(cayleydickson3 a)
{
    return cayleydickson3_make(0, a.data[1], a.data[2], a.data[3],
                               a.data[4], a.data[5], a.data[6], a.data[7]);
}

static inline cayleydickson3
cayleydickson3_from_re(double a)
{
    return cayleydickson3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
cayleydickson3_to_re(cayleydickson3 a)
{
    return a.data[0];
}

static inline cayleydickson3
cayleydickson3_one(void)
{
    return cayleydickson3_from_re(1);
}

static inline cayleydickson3
cayleydickson3_mul(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return octonion_mul(a, b);
    }
    vec2 k01 = {{k.data[0], k.data[1]}};
    cayleydickson2 x0 = {{a.data[0], a.data[1], a.data[2], a.data[3]}};
    cayleydickson2 x1 = {{a.data[4], a.data[5], a.data[6], a.data[7]}};
    cayleydickson2 y0 = {{b.data[0], b.data[1], b.data[2], b.data[3]}};
    cayleydickson2 y1 = {{b.data[4], b.data[5], b.data[6], b.data[7]}};
    cayleydickson2 left = vec4_add(
        cayleydickson2_mul(k01, x0, y0),
        vec4_scale(cayleydickson2_mul(
                       k01, cayleydickson2_conjugate(y1), x1),
                   k.data[2]));
    cayleydickson2 right = vec4_add(
        cayleydickson2_mul(k01, y1, x0),
        cayleydickson2_mul(k01, x1, cayleydickson2_conjugate(y0)));
    return cayleydickson3_make(left.data[0], left.data[1],
                               left.data[2], left.data[3],
                               right.data[0], right.data[1],
                               right.data[2], right.data[3]);
}

static inline cayleydickson3
cayleydickson3_conjugate(cayleydickson3 a)
{
    return vec8_sub(cayleydickson3_re(a), cayleydickson3_im(a));
}

static inline double
cayleydickson3_norm2(vec3 k, cayleydickson3 a)
{
    return a.data[0]*a.data[0]
        - k.data[0]*a.data[1]*a.data[1]
        - k.data[1]*a.data[2]*a.data[2]
        + k.data[0]*k.data[1]*a.data[3]*a.data[3]
        - k.data[2]*a.data[4]*a.data[4]
        + k.data[0]*k.data[2]*a.data[5]*a.data[5]
        + k.data[1]*k.data[2]*a.data[6]*a.data[6]
        - k.data[0]*k.data[1]*k.data[2]*a.data[7]*a.data[7];
}

static inline cayleydickson3
cayleydickson3_inv(vec3 k, cayleydickson3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return octonion_inv(a);
    }
    return vec8_div_scalar(cayleydickson3_conjugate(a),
                           cayleydickson3_norm2(k, a));
}

static inline cayleydickson3
cayleydickson3_div(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    return cayleydickson3_mul(k, a, cayleydickson3_inv(k, b));
}

static inline void
cayleydickson3_basis(cayleydickson3 result[ALGPY_CAYLEYDICKSON3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON3_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson3_coeffs(
    vec3 k,
    cayleydickson3 result[ALGPY_CAYLEYDICKSON3_N *
                          ALGPY_CAYLEYDICKSON3_N])
{
    cayleydickson3 basis[ALGPY_CAYLEYDICKSON3_N];
    int i, j;
    cayleydickson3_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON3_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON3_N + j] =
                cayleydickson3_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
cayleydickson3_minpoly(vec3 k, cayleydickson3 a, double result[3])
{
    result[0] = cayleydickson3_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

static inline void
cayleydickson3_eivals(vec3 k, cayleydickson3 a, complex result[2])
{
    complex s = complex_sqrt_real(-cayleydickson3_norm2(k, cayleydickson3_im(a)));
    result[0] = complex_make(a.data[0] - s.data[0], -s.data[1]);
    result[1] = complex_make(a.data[0] + s.data[0], +s.data[1]);
}

/* Common algebra API. */
static inline cayleydickson3
cayleydickson3_zero(void)
{
    return vec8_zero();
}

static inline int
cayleydickson3_is_real(cayleydickson3 a)
{
    return vec8_is_real(a);
}

static inline cayleydickson3
cayleydickson3_pos(cayleydickson3 a)
{
    return vec8_pos(a);
}

static inline cayleydickson3
cayleydickson3_neg(cayleydickson3 a)
{
    return vec8_neg(a);
}

static inline cayleydickson3
cayleydickson3_add(cayleydickson3 a, cayleydickson3 b)
{
    return vec8_add(a, b);
}

static inline cayleydickson3
cayleydickson3_sub(cayleydickson3 a, cayleydickson3 b)
{
    return vec8_sub(a, b);
}

static inline int
cayleydickson3_equal(cayleydickson3 a, cayleydickson3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
cayleydickson3_not_equal(cayleydickson3 a, cayleydickson3 b)
{
    return !cayleydickson3_equal(a, b);
}

/* Every element satisfies a^2 - 2*re(a)*a + norm2(a) = 0.
 * For a = w + v this gives v^2 = q: complex if q < 0, split if q > 0,
 * and f(a) = f(w) + f'(w)*v if q == 0. */
static inline cayleydickson3
cayleydickson3_slice(cayleydickson3 a, double real, double scale)
{
    cayleydickson3 result;
    int i;
    result.data[0] = real;
    for (i = 1; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : scale*a.data[i];
    }
    return result;
}

static inline cayleydickson3
cayleydickson3_real_slice(vec3 k, complex value)
{
    cayleydickson3 result = cayleydickson3_from_re(value.data[0]);
    int i, p;
    if (value.data[1] == 0) {
        return result;
    }
    for (i = 1; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        double square = -1;
        for (p = 0; p < 3; p++) {
            if (i & (1 << p)) {
                square *= -k.data[p];
            }
        }
        if (square < 0) {
            result.data[i] = value.data[1]/sqrt(-square);
            return result;
        }
    }
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_apply(vec3 k, cayleydickson3 a,
    double (*f)(double), double (*deriv1)(double), complex (*cf)(complex))
{
    double w = a.data[0], q, radius, plus, minus;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return octonion_apply(a, cf);
    }
    if (vec8_is_real(a)) {
        return cayleydickson3_real_slice(k, cf(complex_make(w, 0)));
    }
    q = -cayleydickson3_norm2(k, cayleydickson3_im(a));
    if (q < 0) {
        complex value;
        cayleydickson3 result;
        radius = sqrt(-q);
        value = cf(complex_make(w, radius));
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON3_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = f(w + radius);
        minus = f(w - radius);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
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
                return cayleydickson3_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson3_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        double real = f(w), scale = deriv1(w);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson3_slice(a, real, scale);
        }
    }
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_exp(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_exp, double_exp_deriv1, complex_exp);
}

static inline cayleydickson3
cayleydickson3_log(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_log, double_log_deriv1, complex_log);
}

static inline cayleydickson3
cayleydickson3_sqrt(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_sqrt, double_sqrt_deriv1, complex_sqrt);
}

static inline cayleydickson3
cayleydickson3_cos(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_cos, double_cos_deriv1, complex_cos);
}

static inline cayleydickson3
cayleydickson3_sin(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_sin, double_sin_deriv1, complex_sin);
}

static inline cayleydickson3
cayleydickson3_tan(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_tan, double_tan_deriv1, complex_tan);
}

static inline cayleydickson3
cayleydickson3_acos(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_acos, double_acos_deriv1, complex_acos);
}

static inline cayleydickson3
cayleydickson3_asin(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_asin, double_asin_deriv1, complex_asin);
}

static inline cayleydickson3
cayleydickson3_atan(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_atan, double_atan_deriv1, complex_atan);
}

static inline cayleydickson3
cayleydickson3_cosh(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_cosh, double_cosh_deriv1, complex_cosh);
}

static inline cayleydickson3
cayleydickson3_sinh(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_sinh, double_sinh_deriv1, complex_sinh);
}

static inline cayleydickson3
cayleydickson3_tanh(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_tanh, double_tanh_deriv1, complex_tanh);
}

static inline cayleydickson3
cayleydickson3_acosh(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_acosh, double_acosh_deriv1, complex_acosh);
}

static inline cayleydickson3
cayleydickson3_asinh(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_asinh, double_asinh_deriv1, complex_asinh);
}

static inline cayleydickson3
cayleydickson3_atanh(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_atanh, double_atanh_deriv1, complex_atanh);
}

static inline cayleydickson3
cayleydickson3_square(vec3 k, cayleydickson3 a)
{
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return octonion_square(a);
    }
    return cayleydickson3_slice(a, a.data[0]*a.data[0] - cayleydickson3_norm2(k, cayleydickson3_im(a)), 2*a.data[0]);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline cayleydickson3
cayleydickson3_ipow(vec3 k, cayleydickson3 a, int64_t b)
{
    cayleydickson3 result = cayleydickson3_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson3_mul(k, a, a);
    }
    if (b == -1) {
        return cayleydickson3_inv(k, a);
    }
    if (b < 0) {
        a = cayleydickson3_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = cayleydickson3_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = cayleydickson3_mul(k, a, a);
        }
    }
    return result;
}

static inline cayleydickson3
cayleydickson3_lpows(vec3 k, cayleydickson3 a, double b)
{
    double w = a.data[0], q, radius, plus, minus, real, scale;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1 && k.data[2] == -1) {
        return octonion_lpows(a, b);
    }
    if (b == 0) {
        return cayleydickson3_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson3_square(k, a);
    }
    if (vec8_is_real(a)) {
        real = pow(w, b);
        if (!isnan(real)) {
            return cayleydickson3_from_re(real);
        }
        return cayleydickson3_real_slice(k, complex_lpows(complex_make(w, 0), b));
    }
    q = -cayleydickson3_norm2(k, cayleydickson3_im(a));
    if (q < 0) {
        complex value;
        cayleydickson3 result;
        radius = sqrt(-q);
        value = complex_lpows(complex_make(w, radius), b);
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON3_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = pow(w + radius, b);
        minus = pow(w - radius, b);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
        }
        if (radius < 1e-8*fabs(w)/fmax(1, fabs(b))) {
            real = pow(w, b);
            scale = b*pow(w, b - 1);
            if (isfinite(real) && isfinite(scale)) {
                return cayleydickson3_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson3_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        real = pow(w, b);
        scale = b*pow(w, b - 1);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson3_slice(a, real, scale);
        }
    }
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_rpow(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    if (vec8_is_real(b)) {
        return cayleydickson3_lpows(k, a, b.data[0]);
    }
    return cayleydickson3_exp(k, cayleydickson3_mul(k, cayleydickson3_log(k, a), b));
}

static inline cayleydickson3
cayleydickson3_pow(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    return cayleydickson3_rpow(k, a, b);
}

static inline cayleydickson3
cayleydickson3_lpow(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    return cayleydickson3_exp(k, cayleydickson3_mul(k, b, cayleydickson3_log(k, a)));
}



static inline cayleydickson3
cayleydickson3_radd(double a, cayleydickson3 b)
{
    return cayleydickson3_add(cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_ladd(cayleydickson3 a, double b)
{
    return cayleydickson3_add(a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rsub(double a, cayleydickson3 b)
{
    return cayleydickson3_sub(cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_lsub(cayleydickson3 a, double b)
{
    return cayleydickson3_sub(a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rmul(vec3 k, double a, cayleydickson3 b)
{
    return cayleydickson3_mul(k, cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_lmul(vec3 k, cayleydickson3 a, double b)
{
    return cayleydickson3_mul(k, a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rdiv(vec3 k, double a, cayleydickson3 b)
{
    return cayleydickson3_div(k, cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_ldiv(vec3 k, cayleydickson3 a, double b)
{
    return cayleydickson3_div(k, a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rpows(vec3 k, double a, cayleydickson3 b)
{
    return cayleydickson3_rpow(k, cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_expm1(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_expm1, double_expm1_deriv1, complex_expm1);
}

static inline cayleydickson3
cayleydickson3_log1p(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_apply(k, a, double_log1p, double_log1p_deriv1, complex_log1p);
}

static inline cayleydickson3
cayleydickson3_mix(vec3 k, cayleydickson3 a, cayleydickson3 b, double t)
{
    cayleydickson3 delta = cayleydickson3_mul(k, cayleydickson3_inv(k, a), b);
    return cayleydickson3_mul(k, a, cayleydickson3_lpows(k, delta, t));
}

static inline complex
cayleydickson3_norm(vec3 k, cayleydickson3 a)
{
    return complex_root_real(cayleydickson3_norm2(k, a), 2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
cayleydickson3_lmat(vec3 k, cayleydickson3 a, double result[ALGPY_CAYLEYDICKSON3_N*ALGPY_CAYLEYDICKSON3_N])
{
    cayleydickson3 basis = {{0}};
    cayleydickson3 product;
    int row, col;
    for (col = 0; col < ALGPY_CAYLEYDICKSON3_N; col++) {
        basis.data[col] = 1;
        product = cayleydickson3_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CAYLEYDICKSON3_N; row++) {
            result[row * ALGPY_CAYLEYDICKSON3_N + col] = product.data[row];
        }
    }
}

static inline void
cayleydickson3_rmat(vec3 k, cayleydickson3 a, double result[ALGPY_CAYLEYDICKSON3_N*ALGPY_CAYLEYDICKSON3_N])
{
    cayleydickson3 basis = {{0}};
    cayleydickson3 product;
    int row, col;
    for (row = 0; row < ALGPY_CAYLEYDICKSON3_N; row++) {
        basis.data[row] = 1;
        product = cayleydickson3_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CAYLEYDICKSON3_N; col++) {
            result[row * ALGPY_CAYLEYDICKSON3_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_CAYLEYDICKSON3_H */
