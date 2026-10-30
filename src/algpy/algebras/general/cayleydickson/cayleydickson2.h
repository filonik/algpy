#ifndef ALGPY_CAYLEYDICKSON2_H
#define ALGPY_CAYLEYDICKSON2_H

#define ALGPY_CAYLEYDICKSON2_N 4

#include <stdint.h>
#include <math.h>

#include "../../../fields/complex.h"
#include "../../../fields/double.h"
#include "../../cayleydickson/quaternion.h"
#include "../../../vectorspaces/vec2.h"
#include "../../../vectorspaces/vec4.h"

typedef vec4 cayleydickson2;

static inline cayleydickson2
cayleydickson2_make(double w, double x, double y, double z)
{
    cayleydickson2 r = {{w, x, y, z}};
    return r;
}

static inline cayleydickson2
cayleydickson2_re(cayleydickson2 a)
{
    return cayleydickson2_make(a.data[0], 0, 0, 0);
}

static inline cayleydickson2
cayleydickson2_im(cayleydickson2 a)
{
    return cayleydickson2_make(0, a.data[1], a.data[2], a.data[3]);
}

static inline cayleydickson2
cayleydickson2_from_re(double a)
{
    return cayleydickson2_make(a, 0, 0, 0);
}

static inline double
cayleydickson2_to_re(cayleydickson2 a)
{
    return a.data[0];
}

static inline cayleydickson2
cayleydickson2_one(void)
{
    return cayleydickson2_from_re(1);
}

static inline cayleydickson2
cayleydickson2_mul(vec2 h, cayleydickson2 a, cayleydickson2 b)
{
    if (h.data[0] == -1 && h.data[1] == -1) {
        return quaternion_mul(a, b);
    }
    const double *k = h.data, *x = a.data, *y = b.data;
    return cayleydickson2_make(
      x[0]*y[0] + k[0]*x[1]*y[1] + k[1]*x[2]*y[2] - k[0]*k[1]*x[3]*y[3],
      x[0]*y[1] + x[1]*y[0] - k[1]*x[2]*y[3] + k[1]*x[3]*y[2],
      x[0]*y[2] + k[0]*x[1]*y[3] + x[2]*y[0] - k[0]*x[3]*y[1],
      x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0]);
}

static inline void
cayleydickson2_basis(cayleydickson2 result[ALGPY_CAYLEYDICKSON2_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON2_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson2_coeffs(
    vec2 h,
    cayleydickson2 result[ALGPY_CAYLEYDICKSON2_N * ALGPY_CAYLEYDICKSON2_N])
{
    cayleydickson2 basis[ALGPY_CAYLEYDICKSON2_N];
    int i, j;
    cayleydickson2_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON2_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON2_N + j] =
                cayleydickson2_mul(h, basis[i], basis[j]);
        }
    }
}

static inline cayleydickson2
cayleydickson2_conjugate(cayleydickson2 a)
{
    return vec4_sub(cayleydickson2_re(a), cayleydickson2_im(a));
}

/* Real norm a*conj(a) = w^2 - k0*x^2 - k1*y^2 + k0*k1*z^2 */
static inline double
cayleydickson2_norm2(vec2 k, cayleydickson2 a)
{
    return a.data[0]*a.data[0]
        - k.data[0]*a.data[1]*a.data[1]
        - k.data[1]*a.data[2]*a.data[2]
        + k.data[0]*k.data[1]*a.data[3]*a.data[3];
}

static inline cayleydickson2
cayleydickson2_inv(vec2 k, cayleydickson2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return quaternion_inv(a);
    }
    return vec4_div_scalar(cayleydickson2_conjugate(a), cayleydickson2_norm2(k, a));
}

static inline cayleydickson2
cayleydickson2_div(vec2 k, cayleydickson2 a, cayleydickson2 b)
{
    return cayleydickson2_mul(k, a, cayleydickson2_inv(k, b));
}

/* Every element satisfies a^2 - 2*re(a)*a + norm2(a) = 0.
 * For a = w + v this gives v^2 = q: complex if q < 0, split if q > 0,
 * and f(a) = f(w) + f'(w)*v if q == 0. */
static inline cayleydickson2
cayleydickson2_slice(cayleydickson2 a, double real, double scale)
{
    cayleydickson2 result;
    int i;
    result.data[0] = real;
    for (i = 1; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : scale*a.data[i];
    }
    return result;
}

static inline cayleydickson2
cayleydickson2_real_slice(vec2 k, complex value)
{
    cayleydickson2 result = cayleydickson2_from_re(value.data[0]);
    int i, p;
    if (value.data[1] == 0) {
        return result;
    }
    for (i = 1; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        double square = -1;
        for (p = 0; p < 2; p++) {
            if (i & (1 << p)) {
                square *= -k.data[p];
            }
        }
        if (square < 0) {
            result.data[i] = value.data[1]/sqrt(-square);
            return result;
        }
    }
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_apply(vec2 k, cayleydickson2 a,
    double (*f)(double), double (*deriv1)(double), complex (*cf)(complex))
{
    double w = a.data[0], q, radius, plus, minus;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1) {
        return quaternion_apply(a, cf);
    }
    if (vec4_is_real(a)) {
        return cayleydickson2_real_slice(k, cf(complex_make(w, 0)));
    }
    q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    if (q < 0) {
        complex value;
        cayleydickson2 result;
        radius = sqrt(-q);
        value = cf(complex_make(w, radius));
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON2_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = f(w + radius);
        minus = f(w - radius);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson2_make(NAN, NAN, NAN, NAN);
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
                return cayleydickson2_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson2_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        double real = f(w), scale = deriv1(w);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson2_slice(a, real, scale);
        }
    }
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_exp(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_exp, double_exp_deriv1, complex_exp);
}

static inline cayleydickson2
cayleydickson2_log(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_log, double_log_deriv1, complex_log);
}

static inline cayleydickson2
cayleydickson2_square(vec2 k, cayleydickson2 a)
{
    if (k.data[0] == -1 && k.data[1] == -1) {
        return quaternion_square(a);
    }
    return cayleydickson2_slice(a, a.data[0]*a.data[0] - cayleydickson2_norm2(k, cayleydickson2_im(a)), 2*a.data[0]);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline cayleydickson2
cayleydickson2_ipow(vec2 k, cayleydickson2 a, int64_t b)
{
    cayleydickson2 result = cayleydickson2_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson2_mul(k, a, a);
    }
    if (b == -1) {
        return cayleydickson2_inv(k, a);
    }
    if (b < 0) {
        a = cayleydickson2_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = cayleydickson2_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = cayleydickson2_mul(k, a, a);
        }
    }
    return result;
}

static inline cayleydickson2
cayleydickson2_lpows(vec2 k, cayleydickson2 a, double b)
{
    double w = a.data[0], q, radius, plus, minus, real, scale;
    int i;
    if (k.data[0] == -1 && k.data[1] == -1) {
        return quaternion_lpows(a, b);
    }
    if (b == 0) {
        return cayleydickson2_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cayleydickson2_square(k, a);
    }
    if (vec4_is_real(a)) {
        real = pow(w, b);
        if (!isnan(real)) {
            return cayleydickson2_from_re(real);
        }
        return cayleydickson2_real_slice(k, complex_lpows(complex_make(w, 0), b));
    }
    q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    if (q < 0) {
        complex value;
        cayleydickson2 result;
        radius = sqrt(-q);
        value = complex_lpows(complex_make(w, radius), b);
        result.data[0] = value.data[0];
        for (i = 1; i < ALGPY_CAYLEYDICKSON2_N; i++) {
            result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
        }
        return result;
    }
    if (q > 0) {
        radius = sqrt(q);
        plus = pow(w + radius, b);
        minus = pow(w - radius, b);
        if (isnan(plus) || isnan(minus)) {
            return cayleydickson2_make(NAN, NAN, NAN, NAN);
        }
        if (radius < 1e-8*fabs(w)/fmax(1, fabs(b))) {
            real = pow(w, b);
            scale = b*pow(w, b - 1);
            if (isfinite(real) && isfinite(scale)) {
                return cayleydickson2_slice(a, real, scale);
            }
        }
        if (w + radius != w - radius) {
            return cayleydickson2_slice(a, 0.5*plus + 0.5*minus, (0.5*plus - 0.5*minus)/radius);
        }
    }
    if (q == 0 || q > 0) {
        real = pow(w, b);
        scale = b*pow(w, b - 1);
        if (isfinite(real) && isfinite(scale)) {
            return cayleydickson2_slice(a, real, scale);
        }
    }
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_rpow(vec2 k, cayleydickson2 a, cayleydickson2 b)
{
    if (vec4_is_real(b)) {
        return cayleydickson2_lpows(k, a, b.data[0]);
    }
    return cayleydickson2_exp(k, cayleydickson2_mul(k, cayleydickson2_log(k, a), b));
}

static inline cayleydickson2
cayleydickson2_pow(vec2 k, cayleydickson2 a, cayleydickson2 b)
{
    return cayleydickson2_rpow(k, a, b);
}

static inline cayleydickson2
cayleydickson2_lpow(vec2 k, cayleydickson2 a, cayleydickson2 b)
{
    return cayleydickson2_exp(k, cayleydickson2_mul(k, b, cayleydickson2_log(k, a)));
}

static inline cayleydickson2
cayleydickson2_rpows(vec2 k, double a, cayleydickson2 b)
{
    return cayleydickson2_rpow(k, cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_mix(vec2 k, cayleydickson2 a, cayleydickson2 b, double t)
{
    cayleydickson2 delta = cayleydickson2_mul(k, cayleydickson2_inv(k, a), b);
    return cayleydickson2_mul(k, a, cayleydickson2_lpows(k, delta, t));
}

static inline cayleydickson2
cayleydickson2_sqrt(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_sqrt, double_sqrt_deriv1, complex_sqrt);
}

static inline cayleydickson2
cayleydickson2_cos(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_cos, double_cos_deriv1, complex_cos);
}

static inline cayleydickson2
cayleydickson2_sin(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_sin, double_sin_deriv1, complex_sin);
}

static inline cayleydickson2
cayleydickson2_tan(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_tan, double_tan_deriv1, complex_tan);
}

static inline cayleydickson2
cayleydickson2_cosh(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_cosh, double_cosh_deriv1, complex_cosh);
}

static inline cayleydickson2
cayleydickson2_sinh(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_sinh, double_sinh_deriv1, complex_sinh);
}

static inline cayleydickson2
cayleydickson2_tanh(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_tanh, double_tanh_deriv1, complex_tanh);
}

/* Minimal polynomial t^2 - 2*w*t + norm2(a) (coefficients in ascending order) */
static inline void
cayleydickson2_minpoly(vec2 k, cayleydickson2 a, double result[3])
{
    result[0] = cayleydickson2_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/* Roots w +/- sqrt(-norm2(im(a))) of the minimal polynomial, as complex numbers. */
static inline void
cayleydickson2_eivals(vec2 k, cayleydickson2 a, complex result[2])
{
    complex s = complex_sqrt_real(-cayleydickson2_norm2(k, cayleydickson2_im(a)));
    result[0] = complex_make(a.data[0] - s.data[0], -s.data[1]);
    result[1] = complex_make(a.data[0] + s.data[0], +s.data[1]);
}

/* Common algebra API. */
static inline cayleydickson2
cayleydickson2_zero(void)
{
    return vec4_zero();
}

static inline int
cayleydickson2_is_real(cayleydickson2 a)
{
    return vec4_is_real(a);
}

static inline cayleydickson2
cayleydickson2_pos(cayleydickson2 a)
{
    return vec4_pos(a);
}

static inline cayleydickson2
cayleydickson2_neg(cayleydickson2 a)
{
    return vec4_neg(a);
}

static inline cayleydickson2
cayleydickson2_add(cayleydickson2 a, cayleydickson2 b)
{
    return vec4_add(a, b);
}

static inline cayleydickson2
cayleydickson2_sub(cayleydickson2 a, cayleydickson2 b)
{
    return vec4_sub(a, b);
}

static inline int
cayleydickson2_equal(cayleydickson2 a, cayleydickson2 b)
{
    return vec4_all_eq(a, b);
}

static inline int
cayleydickson2_not_equal(cayleydickson2 a, cayleydickson2 b)
{
    return !cayleydickson2_equal(a, b);
}

static inline cayleydickson2
cayleydickson2_acos(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_acos, double_acos_deriv1, complex_acos);
}

static inline cayleydickson2
cayleydickson2_asin(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_asin, double_asin_deriv1, complex_asin);
}

static inline cayleydickson2
cayleydickson2_atan(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_atan, double_atan_deriv1, complex_atan);
}

static inline cayleydickson2
cayleydickson2_acosh(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_acosh, double_acosh_deriv1, complex_acosh);
}

static inline cayleydickson2
cayleydickson2_asinh(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_asinh, double_asinh_deriv1, complex_asinh);
}

static inline cayleydickson2
cayleydickson2_atanh(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_atanh, double_atanh_deriv1, complex_atanh);
}



static inline cayleydickson2
cayleydickson2_radd(double a, cayleydickson2 b)
{
    return cayleydickson2_add(cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_ladd(cayleydickson2 a, double b)
{
    return cayleydickson2_add(a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_rsub(double a, cayleydickson2 b)
{
    return cayleydickson2_sub(cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_lsub(cayleydickson2 a, double b)
{
    return cayleydickson2_sub(a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_rmul(vec2 k, double a, cayleydickson2 b)
{
    return cayleydickson2_mul(k, cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_lmul(vec2 k, cayleydickson2 a, double b)
{
    return cayleydickson2_mul(k, a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_rdiv(vec2 k, double a, cayleydickson2 b)
{
    return cayleydickson2_div(k, cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_ldiv(vec2 k, cayleydickson2 a, double b)
{
    return cayleydickson2_div(k, a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_expm1(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_expm1, double_expm1_deriv1, complex_expm1);
}

static inline cayleydickson2
cayleydickson2_log1p(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_apply(k, a, double_log1p, double_log1p_deriv1, complex_log1p);
}

static inline complex
cayleydickson2_norm(vec2 k, cayleydickson2 a)
{
    return complex_root_real(cayleydickson2_norm2(k, a), 2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
cayleydickson2_lmat(vec2 k, cayleydickson2 a, double result[ALGPY_CAYLEYDICKSON2_N*ALGPY_CAYLEYDICKSON2_N])
{
    cayleydickson2 basis = {{0}};
    cayleydickson2 product;
    int row, col;
    for (col = 0; col < ALGPY_CAYLEYDICKSON2_N; col++) {
        basis.data[col] = 1;
        product = cayleydickson2_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CAYLEYDICKSON2_N; row++) {
            result[row * ALGPY_CAYLEYDICKSON2_N + col] = product.data[row];
        }
    }
}

static inline void
cayleydickson2_rmat(vec2 k, cayleydickson2 a, double result[ALGPY_CAYLEYDICKSON2_N*ALGPY_CAYLEYDICKSON2_N])
{
    cayleydickson2 basis = {{0}};
    cayleydickson2 product;
    int row, col;
    for (row = 0; row < ALGPY_CAYLEYDICKSON2_N; row++) {
        basis.data[row] = 1;
        product = cayleydickson2_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CAYLEYDICKSON2_N; col++) {
            result[row * ALGPY_CAYLEYDICKSON2_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_CAYLEYDICKSON2_H */
