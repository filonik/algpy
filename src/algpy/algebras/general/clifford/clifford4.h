#ifndef ALGPY_CLIFFORD4_H
#define ALGPY_CLIFFORD4_H

#define ALGPY_CLIFFORD4_N 16

#include <stdint.h>
#include <math.h>

#include "../../../fields/complex.h"
#include "../../../polynomials/hermite.h"
#include "clifford3.h"
#include "../../../vectorspaces/vec4.h"
#include "../../../vectorspaces/vec16.h"

typedef vec16 clifford4;

static inline int popcount(unsigned int x)
{
#if defined(__clang__) || defined(__GNUC__)
    return __builtin_popcount(x);
#else
    int count = 0;
    while (x != 0) {
        x &= x - 1;  /* Clear the lowest set bit. */
        ++count;
    }
    return count;
#endif
}

static inline clifford4
clifford4_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012,
                    double x3, double x03, double x13, double x013,
                    double x23, double x023, double x123, double x0123)
{
    clifford4 r = {{x, x0, x1, x01, x2, x02, x12, x012,
                         x3, x03, x13, x013, x23, x023, x123, x0123}};
    return r;
}

static inline clifford4
clifford4_grade0(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline clifford4
clifford4_grade1(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0);
}

static inline clifford4
clifford4_grade2(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, 0, 0, x[3], 0, x[5], x[6], 0,
                          0, x[9], x[10], 0, x[12], 0, 0, 0);
}

static inline clifford4
clifford4_grade3(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, 0, 0, 0, 0, 0, 0, x[7],
                          0, 0, 0, x[11], 0, x[13], x[14], 0);
}

static inline clifford4
clifford4_grade4(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, x[15]);
}

static inline clifford4
clifford4_from_re(double a)
{
    return clifford4_make(a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
clifford4_to_re(clifford4 a)
{
    return a.data[0];
}

static inline clifford4
clifford4_one(void)
{
    return clifford4_from_re(1);
}

static inline clifford4
clifford4_radd(double a, clifford4 b) {
  return vec16_add(clifford4_from_re(a), b);
}

static inline clifford4
clifford4_ladd(clifford4 a, double b) {
  return vec16_add(a, clifford4_from_re(b));
}

static inline clifford4
clifford4_rsub(double a, clifford4 b) {
  return vec16_sub(clifford4_from_re(a), b);
}

static inline clifford4
clifford4_lsub(clifford4 a, double b) {
  return vec16_sub(a, clifford4_from_re(b));
}

static inline clifford4
clifford4_mul(vec4 k, clifford4 a, clifford4 b)
{
    clifford4 r = {{0}};
    int i, j, p, q;
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD4_N; j++) {
            double factor = 1;
            int swaps = 0;
            for (p = 0; p < 4; p++) {
                if (!(j & (1 << p))) {
                    continue;
                }
                for (q = p + 1; q < 4; q++) {
                    if (i & (1 << q)) {
                        swaps++;
                    }
                }
                if (i & (1 << p)) {
                    factor *= k.data[p];
                }
            }
            if (swaps & 1) {
                factor = -factor;
            }
            r.data[i ^ j] += factor * a.data[i] * b.data[j];
        }
    }
    return r;
}

/* Clifford conjugate: grade g is scaled by (-1)^(g(g+1)/2). */
static inline clifford4
clifford4_conjugate(clifford4 a)
{
    clifford4 r = a;
    int i;
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        int g = popcount(i);
        if (((g * (g + 1)) / 2) & 1) {
            r.data[i] = -r.data[i];
        }
    }
    return r;
}

static inline void
clifford4_basis(clifford4 result[ALGPY_CLIFFORD4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        double components[ALGPY_CLIFFORD4_N] = {0};
        components[i] = 1;
        for (j = 0; j < ALGPY_CLIFFORD4_N; j++) {
            result[i].data[j] = components[j];
        }
    }
}

static inline void
clifford4_coeffs(
    vec4 k,
    clifford4 result[ALGPY_CLIFFORD4_N * ALGPY_CLIFFORD4_N])
{
    clifford4 basis[ALGPY_CLIFFORD4_N];
    int i, j;
    clifford4_basis(basis);
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD4_N; j++) {
            result[i * ALGPY_CLIFFORD4_N + j] =
                clifford4_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
clifford4_minpoly(vec4 k, clifford4 a, double result[5])
{
    const int n = 4;
    clifford4 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * clifford4_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            clifford4 adj = clifford4_lsub(uk, ck);
            uk = clifford4_mul(k, a, adj);
        }
    }
}

/* Norm (degree 4), the constant coefficient of the generic minimal polynomial. */
static inline double
clifford4_norm4(vec4 k, clifford4 a)
{
    double p[5];
    clifford4_minpoly(k, a, p);
    return p[0];
}

static inline complex
clifford4_norm(vec4 k, clifford4 a)
{
    return complex_root_real(clifford4_norm4(k, a), 4);
}

/* Inverse from the Cayley-Hamilton theorem: p(a) = 0 with p(0) = norm4. */
static inline clifford4
clifford4_inv(vec4 k, clifford4 a)
{
    double p[5];
    clifford4 a2, a3, r;
    clifford4_minpoly(k, a, p);
    a2 = clifford4_mul(k, a, a);
    a3 = clifford4_mul(k, a2, a);
    r = vec16_add(vec16_add(a3, vec16_scale(a2, p[3])),
                  vec16_add(vec16_scale(a, p[2]), clifford4_from_re(p[1])));
    return vec16_div_scalar(vec16_neg(r), p[0]);
}

static inline clifford4
clifford4_div(vec4 k, clifford4 a, clifford4 b)
{
    return clifford4_mul(k, a, clifford4_inv(k, b));
}

static inline void
clifford4_eivals(vec4 k, clifford4 a, complex result[4])
{
    double p[5];
    clifford4_minpoly(k, a, p);
    poly4_coeffs_to_roots(p, result);
}

/* Common algebra API. */
static inline clifford4
clifford4_zero(void)
{
    return vec16_zero();
}

static inline int
clifford4_is_real(clifford4 a)
{
    return vec16_is_real(a);
}

static inline clifford4
clifford4_pos(clifford4 a)
{
    return vec16_pos(a);
}

static inline clifford4
clifford4_neg(clifford4 a)
{
    return vec16_neg(a);
}

static inline clifford4
clifford4_add(clifford4 a, clifford4 b)
{
    return vec16_add(a, b);
}

static inline clifford4
clifford4_sub(clifford4 a, clifford4 b)
{
    return vec16_sub(a, b);
}

static inline int
clifford4_equal(clifford4 a, clifford4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
clifford4_not_equal(clifford4 a, clifford4 b)
{
    return !clifford4_equal(a, b);
}

static inline clifford4
clifford4_re(clifford4 a)
{
    return clifford4_from_re(a.data[0]);
}

static inline clifford4
clifford4_im(clifford4 a)
{
    return clifford4_sub(a, clifford4_re(a));
}

/* Scalar part of a*b: only equal basis blades contribute. */
static inline double
clifford4_scalar_product(vec4 k, clifford4 a, clifford4 b)
{
    const double *x = a.data, *y = b.data;
    double k0 = k.data[0], k1 = k.data[1], k2 = k.data[2], k3 = k.data[3];
    return x[0]*y[0] + k0*x[1]*y[1] + k1*x[2]*y[2] - k0*k1*x[3]*y[3]
        + k2*x[4]*y[4] - k0*k2*x[5]*y[5] - k1*k2*x[6]*y[6] - k0*k1*k2*x[7]*y[7]
        + k3*x[8]*y[8] - k0*k3*x[9]*y[9] - k1*k3*x[10]*y[10] - k0*k1*k3*x[11]*y[11]
        - k2*k3*x[12]*y[12] - k0*k2*k3*x[13]*y[13] - k1*k2*k3*x[14]*y[14]
        + k0*k1*k2*k3*x[15]*y[15];
}

/* Principal complex continuation for scalar inputs, in the first available
 * negative-square basis blade, as in the other Clifford orders. */
static inline clifford4
clifford4_real_slice(vec4 k, complex value)
{
    double k0 = k.data[0], k1 = k.data[1], k2 = k.data[2], k3 = k.data[3];
    double squares[16] = {1, k0, k1, -k0*k1, k2, -k0*k2, -k1*k2, -k0*k1*k2,
        k3, -k0*k3, -k1*k3, -k0*k1*k3, -k2*k3, -k0*k2*k3, -k1*k2*k3, k0*k1*k2*k3};
    clifford4 result = clifford4_from_re(value.data[0]);
    int i;
    if (value.data[1] == 0) {
        return result;
    }
    for (i = 1; i < 16; i++) {
        if (squares[i] < 0) {
            result.data[i] = value.data[1]/sqrt(-squares[i]);
            return result;
        }
    }
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                          NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* Any unused generator leaves a three-generator Clifford subalgebra. */
static inline int
clifford4_unused_generator(clifford4 a)
{
    int generator, i;
    for (generator = 0; generator < 4; generator++) {
        int used = 0;
        for (i = 1; i < 16; i++) {
            if ((i & (1 << generator)) && a.data[i] != 0) {
                used = 1;
                break;
            }
        }
        if (!used) {
            return generator;
        }
    }
    return -1;
}

static inline clifford4
clifford4_apply_subalgebra(vec4 k, clifford4 a, int generator,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    vec3 subk;
    clifford3 suba, value;
    clifford4 result = {{0}};
    int low = (1 << generator) - 1, i, j = 0;
    for (i = 0; i < 4; i++) {
        if (i != generator) {
            subk.data[j++] = k.data[i];
        }
    }
    for (i = 0; i < 8; i++) {
        int mask = (i & low) | ((i & ~low) << 1);
        suba.data[i] = a.data[mask];
    }
    value = clifford3_apply(subk, suba, f, deriv1, deriv2, deriv3, cf, cderiv1);
    if (isnan(value.data[0])) {
        return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                              NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
    }
    for (i = 0; i < 8; i++) {
        int mask = (i & low) | ((i & ~low) << 1);
        result.data[mask] = value.data[i];
    }
    return result;
}

/* Find b^2 = constant + linear*b, allowing only multiplication roundoff.
 * The scalar coefficient of b is zero. This also recognizes idempotents
 * whose two eigenvalues have unequal generic multiplicities. */
static inline int
clifford4_quadratic_relation(clifford4 b, clifford4 b2, double *linear)
{
    int i, pivot = 1;
    double bound, size = 0;
    for (i = 2; i < 16; i++) {
        if (fabs(b.data[i]) > fabs(b.data[pivot])) {
            pivot = i;
        }
    }
    *linear = b2.data[pivot]/b.data[pivot];
    if (!isfinite(*linear)) {
        return 0;
    }
    for (i = 0; i < 16; i++) {
        size = fmax(size, fabs(b2.data[i]));
    }
    bound = 32*DBL_EPSILON*(size + fabs(*linear*b.data[pivot]));
    for (i = 1; i < 16; i++) {
        double value = *linear*b.data[i];
        if (fabs(b2.data[i] - value) > bound) {
            return 0;
        }
    }
    return 1;
}

static inline clifford4
clifford4_apply_quadratic(double w, clifford4 b, double constant, double linear,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    double shift = linear/2;
    if (hermite_quadratic_apply(w + shift, constant + shift*shift, 0,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        clifford4 result = vec16_scale(b, value.odd);
        result.data[0] = value.even - shift*value.odd;
        return result;
    }
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                          NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* Detect a cubic element relation before imposing unnecessary derivatives
 * from the generic quartic. Two coefficient pivots determine the relation;
 * every remaining coefficient must agree to multiplication roundoff. */
static inline int
clifford4_cubic_relation(clifford4 b, clifford4 b2, clifford4 b3, double p[4])
{
    int i, first = 1, second = 1;
    double ratio, largest = 0, linear, quadratic, bound = 0;
    for (i = 2; i < 16; i++) {
        if (fabs(b.data[i]) > fabs(b.data[first])) {
            first = i;
        }
    }
    ratio = b2.data[first]/b.data[first];
    for (i = 1; i < 16; i++) {
        double residual = fabs(b2.data[i] - ratio*b.data[i]);
        if (residual > largest) {
            largest = residual;
            second = i;
        }
    }
    if (largest == 0) {
        return 0;
    }
    quadratic = (b3.data[second] - b3.data[first]/b.data[first]*b.data[second])/
        (b2.data[second] - ratio*b.data[second]);
    linear = (b3.data[first] - quadratic*b2.data[first])/b.data[first];
    if (!isfinite(linear) || !isfinite(quadratic)) {
        return 0;
    }
    for (i = 0; i < 16; i++) {
        bound = fmax(bound, fabs(b3.data[i]) + fabs(quadratic*b2.data[i]) + fabs(linear*b.data[i]));
    }
    bound *= 64*DBL_EPSILON;
    for (i = 1; i < 16; i++) {
        double q = quadratic*b2.data[i], l = linear*b.data[i];
        if (fabs(b3.data[i] - q - l) > bound) {
            return 0;
        }
    }
    p[0] = quadratic*b2.data[0] - b3.data[0];
    p[1] = -linear;
    p[2] = -quadratic;
    p[3] = 1;
    return 1;
}

/* For b=a-w, s_j=<b^j>_0, the centred generic minimal polynomial is
 * b^4 - 2*s2*b^2 - (4/3)*s3*b + 2*s2^2-s4 = 0.
 * Only b^2 and b^3 require full products: s4=<b^2*b^2>_0 is a contraction.
 * Scalar, subalgebra and lower-degree element relations are separate paths.
 * The polynomial engine knows only the quartic and the field callbacks. */
static inline clifford4
clifford4_apply(vec4 k, clifford4 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double w = a.data[0], linear, coefficients[4], cubic[4];
    clifford4 b = a, b2, b3, result;
    int generator, ok, i;
    if (vec16_is_real(a)) {
        return clifford4_real_slice(k, cf(complex_make(w, 0)));
    }
    generator = clifford4_unused_generator(a);
    if (generator >= 0) {
        return clifford4_apply_subalgebra(k, a, generator,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    b.data[0] = 0;
    b2 = clifford4_mul(k, b, b);
    if (clifford4_quadratic_relation(b, b2, &linear)) {
        return clifford4_apply_quadratic(w, b, b2.data[0], linear,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    b3 = clifford4_mul(k, b2, b);
    if (clifford4_cubic_relation(b, b2, b3, cubic)) {
        complex roots[4];
        poly3_coeffs_to_roots(cubic, roots);
        ok = hermite_interpolate(3, w, roots, f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients);
    } else {
        double s2 = b2.data[0], s3 = b3.data[0];
        double s4 = clifford4_scalar_product(k, b2, b2);
        double quartic[5] = {2*s2*s2 - s4, -4*s3/3, -2*s2, 0, 1};
        ok = poly4_apply(quartic, w, f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients);
    }
    if (!ok) {
        return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                              NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
    }
    for (i = 0; i < 16; i++) {
        result.data[i] = coefficients[1]*b.data[i] + coefficients[2]*b2.data[i] + coefficients[3]*b3.data[i];
    }
    result.data[0] += coefficients[0];
    return result;
}

static inline clifford4
clifford4_exp(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_exp, double_exp_deriv1,
        double_exp_deriv2, double_exp_deriv3, complex_exp, complex_exp_deriv1);
}

static inline clifford4
clifford4_log(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_log, double_log_deriv1,
        double_log_deriv2, double_log_deriv3, complex_log, complex_log_deriv1);
}

static inline clifford4
clifford4_sqrt(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_sqrt, double_sqrt_deriv1,
        double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt, complex_sqrt_deriv1);
}

static inline clifford4
clifford4_cos(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_cos, double_cos_deriv1,
        double_cos_deriv2, double_cos_deriv3, complex_cos, complex_cos_deriv1);
}

static inline clifford4
clifford4_sin(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_sin, double_sin_deriv1,
        double_sin_deriv2, double_sin_deriv3, complex_sin, complex_sin_deriv1);
}

static inline clifford4
clifford4_tan(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_tan, double_tan_deriv1,
        double_tan_deriv2, double_tan_deriv3, complex_tan, complex_tan_deriv1);
}

static inline clifford4
clifford4_acos(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_acos, double_acos_deriv1,
        double_acos_deriv2, double_acos_deriv3, complex_acos, complex_acos_deriv1);
}

static inline clifford4
clifford4_asin(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_asin, double_asin_deriv1,
        double_asin_deriv2, double_asin_deriv3, complex_asin, complex_asin_deriv1);
}

static inline clifford4
clifford4_atan(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_atan, double_atan_deriv1,
        double_atan_deriv2, double_atan_deriv3, complex_atan, complex_atan_deriv1);
}

static inline clifford4
clifford4_cosh(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_cosh, double_cosh_deriv1,
        double_cosh_deriv2, double_cosh_deriv3, complex_cosh, complex_cosh_deriv1);
}

static inline clifford4
clifford4_sinh(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_sinh, double_sinh_deriv1,
        double_sinh_deriv2, double_sinh_deriv3, complex_sinh, complex_sinh_deriv1);
}

static inline clifford4
clifford4_tanh(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_tanh, double_tanh_deriv1,
        double_tanh_deriv2, double_tanh_deriv3, complex_tanh, complex_tanh_deriv1);
}

static inline clifford4
clifford4_acosh(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_acosh, double_acosh_deriv1,
        double_acosh_deriv2, double_acosh_deriv3, complex_acosh, complex_acosh_deriv1);
}

static inline clifford4
clifford4_asinh(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_asinh, double_asinh_deriv1,
        double_asinh_deriv2, double_asinh_deriv3, complex_asinh, complex_asinh_deriv1);
}

static inline clifford4
clifford4_atanh(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_atanh, double_atanh_deriv1,
        double_atanh_deriv2, double_atanh_deriv3, complex_atanh, complex_atanh_deriv1);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline clifford4
clifford4_ipow(vec4 k, clifford4 a, int64_t b)
{
    clifford4 result = clifford4_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return clifford4_mul(k, a, a);
    }
    if (b == -1) {
        return clifford4_inv(k, a);
    }
    if (b < 0) {
        a = clifford4_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = clifford4_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = clifford4_mul(k, a, a);
        }
    }
    return result;
}

static inline clifford4
clifford4_rpow(vec4 k, clifford4 a, clifford4 b)
{
    return clifford4_exp(k, clifford4_mul(k, clifford4_log(k, a), b));
}

static inline clifford4
clifford4_pow(vec4 k, clifford4 a, clifford4 b)
{
    return clifford4_rpow(k, a, b);
}

static inline clifford4
clifford4_lpow(vec4 k, clifford4 a, clifford4 b)
{
    return clifford4_exp(k, clifford4_mul(k, b, clifford4_log(k, a)));
}

static inline clifford4
clifford4_square(vec4 k, clifford4 a)
{
    return clifford4_mul(k, a, a);
}

static inline clifford4
clifford4_rmul(vec4 k, double a, clifford4 b)
{
    return clifford4_mul(k, clifford4_from_re(a), b);
}

static inline clifford4
clifford4_lmul(vec4 k, clifford4 a, double b)
{
    return clifford4_mul(k, a, clifford4_from_re(b));
}

static inline clifford4
clifford4_rdiv(vec4 k, double a, clifford4 b)
{
    return clifford4_div(k, clifford4_from_re(a), b);
}

static inline clifford4
clifford4_ldiv(vec4 k, clifford4 a, double b)
{
    return clifford4_div(k, a, clifford4_from_re(b));
}

static inline clifford4
clifford4_rpows(vec4 k, double a, clifford4 b)
{
    return clifford4_rpow(k, clifford4_from_re(a), b);
}

static inline clifford4
clifford4_lpows(vec4 k, clifford4 a, double b)
{
    return clifford4_rpow(k, a, clifford4_from_re(b));
}

static inline clifford4
clifford4_expm1(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_expm1, double_expm1_deriv1,
        double_expm1_deriv2, double_expm1_deriv3, complex_expm1, complex_expm1_deriv1);
}

static inline clifford4
clifford4_log1p(vec4 k, clifford4 a)
{
    return clifford4_apply(k, a, double_log1p, double_log1p_deriv1,
        double_log1p_deriv2, double_log1p_deriv3, complex_log1p, complex_log1p_deriv1);
}

static inline clifford4
clifford4_mix(vec4 k, clifford4 a, clifford4 b, double t)
{
    clifford4 delta = clifford4_mul(k, clifford4_inv(k, a), b);
    return clifford4_mul(k, a, clifford4_lpows(k, delta, t));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
clifford4_lmat(vec4 k, clifford4 a, double result[ALGPY_CLIFFORD4_N*ALGPY_CLIFFORD4_N])
{
    clifford4 basis = {{0}};
    clifford4 product;
    int row, col;
    for (col = 0; col < ALGPY_CLIFFORD4_N; col++) {
        basis.data[col] = 1;
        product = clifford4_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CLIFFORD4_N; row++) {
            result[row * ALGPY_CLIFFORD4_N + col] = product.data[row];
        }
    }
}

static inline void
clifford4_rmat(vec4 k, clifford4 a, double result[ALGPY_CLIFFORD4_N*ALGPY_CLIFFORD4_N])
{
    clifford4 basis = {{0}};
    clifford4 product;
    int row, col;
    for (row = 0; row < ALGPY_CLIFFORD4_N; row++) {
        basis.data[row] = 1;
        product = clifford4_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CLIFFORD4_N; col++) {
            result[row * ALGPY_CLIFFORD4_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_CLIFFORD4_H */
