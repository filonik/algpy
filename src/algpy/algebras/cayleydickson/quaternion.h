/*
 * Core quaternion arithmetic with fixed generator squares -1.
 *
 * Pure C with no dependency on Python or NumPy, so that it can be reused and
 * tested independently of the wrapper code in `_quaternion.c`.
 */
#ifndef ALGPY_QUATERNION_H
#define ALGPY_QUATERNION_H

#include <stdint.h>
#include "../../fields/complex.h"
#include "../../vectorspaces/vec2.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_QUATERNION_N ALGPY_VEC4_N

typedef vec4 quaternion;

static inline quaternion
quaternion_make(double x0, double x1, double x2, double x3)
{
    quaternion result = {{x0, x1, x2, x3}};
    return result;
}

static inline quaternion
quaternion_re(quaternion a)
{
    return quaternion_make(a.data[0], 0, 0, 0);
}

static inline quaternion
quaternion_im(quaternion a)
{
    return quaternion_make(0, a.data[1], a.data[2], a.data[3]);
}

static inline quaternion
quaternion_from_re(double a)
{
    return quaternion_make(a, 0, 0, 0);
}

static inline double
quaternion_to_re(quaternion a)
{
    return a.data[0];
}

static inline quaternion
quaternion_one(void)
{
    return quaternion_make(1, 0, 0, 0);
}

/* Hamilton product */
static inline quaternion
quaternion_mul(quaternion a, quaternion b)
{
    if (vec4_is_real(b)) {
        return vec4_scale(a, b.data[0]);
    }
    if (vec4_is_real(a)) {
        return vec4_scale(b, a.data[0]);
    }
    return quaternion_make(
        a.data[0] * b.data[0] - a.data[1] * b.data[1] -
            a.data[2] * b.data[2] - a.data[3] * b.data[3],
        a.data[0] * b.data[1] + a.data[1] * b.data[0] +
            a.data[2] * b.data[3] - a.data[3] * b.data[2],
        a.data[0] * b.data[2] - a.data[1] * b.data[3] +
            a.data[2] * b.data[0] + a.data[3] * b.data[1],
        a.data[0] * b.data[3] + a.data[1] * b.data[2] -
            a.data[2] * b.data[1] + a.data[3] * b.data[0]);
}

static inline void
quaternion_basis(quaternion result[ALGPY_QUATERNION_N])
{
    int i, j;
    for (i = 0; i < ALGPY_QUATERNION_N; i++) {
            for (j = 0; j < ALGPY_QUATERNION_N; j++) {
                result[i].data[j] = i == j ? 1 : 0;
            }
    }
}

static inline void
quaternion_coeffs(quaternion result[ALGPY_QUATERNION_N * ALGPY_QUATERNION_N])
{
    quaternion basis[ALGPY_QUATERNION_N];
    int i, j;
    quaternion_basis(basis);
    for (i = 0; i < ALGPY_QUATERNION_N; i++) {
            for (j = 0; j < ALGPY_QUATERNION_N; j++) {
                result[i * ALGPY_QUATERNION_N + j] =
                    quaternion_mul(basis[i], basis[j]);
            }
    }
}

static inline double
quaternion_norm2(quaternion a)
{
    return a.data[0]*a.data[0] + a.data[1]*a.data[1] + a.data[2]*a.data[2] + a.data[3]*a.data[3];
}

static inline double
quaternion_norm(quaternion a)
{
    return sqrt(quaternion_norm2(a));
}

static inline double
quaternion_vector_norm2(quaternion a)
{
    return a.data[1]*a.data[1] + a.data[2]*a.data[2] + a.data[3]*a.data[3];
}

static inline double
quaternion_vector_norm(quaternion a)
{
    double r0 = hypot(a.data[1], a.data[2]);
    return hypot(r0, a.data[3]);
}

static inline quaternion
quaternion_conjugate(quaternion a)
{
    return quaternion_make(a.data[0], -a.data[1], -a.data[2], -a.data[3]);
}

static inline quaternion
quaternion_inv(quaternion a)
{
    return vec4_div_scalar(quaternion_conjugate(a), quaternion_norm2(a));
}

static inline quaternion
quaternion_div(quaternion a, quaternion b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return quaternion_mul(a, quaternion_inv(b));
}

/* Lift a complex value along the input's imaginary direction. Real inputs
 * use the first imaginary basis vector to select the principal branch. */
static inline quaternion
quaternion_complex_slice(quaternion a, complex value, double radius)
{
    quaternion result = {{0}};
    int i;
    result.data[0] = value.data[0];
    if (radius == 0) {
        result.data[1] = value.data[1];
        return result;
    }
    for (i = 1; i < ALGPY_QUATERNION_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
    }
    return result;
}

static inline quaternion
quaternion_apply(quaternion a, complex (*f)(complex))
{
    double radius = quaternion_vector_norm(a);
    complex value = f(complex_make(a.data[0], radius));
    return quaternion_complex_slice(a, value, radius);
}

static inline quaternion
quaternion_exp(quaternion a)
{
    return quaternion_apply(a, complex_exp);
}

/* Principal logarithm; for negative reals the pure part points along i. */
static inline quaternion
quaternion_log(quaternion a)
{
    return quaternion_apply(a, complex_log);
}

static inline quaternion
quaternion_square(quaternion a)
{
    return quaternion_make(
        a.data[0]*a.data[0] - quaternion_vector_norm2(a), 2*a.data[0]*a.data[1], 2*a.data[0]*a.data[2], 2*a.data[0]*a.data[3]);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline quaternion
quaternion_ipow(quaternion a, int64_t b)
{
    quaternion result = quaternion_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return quaternion_mul(a, a);
    }
    if (b == -1) {
        return quaternion_inv(a);
    }
    if (b < 0) {
        a = quaternion_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = quaternion_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = quaternion_mul(a, a);
        }
    }
    return result;
}

static inline quaternion
quaternion_lpows(quaternion a, double b)
{
    double radius, real;
    complex value;
    if (b == 0) {
        return quaternion_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return quaternion_square(a);
    }
    radius = quaternion_vector_norm(a);
    if (radius == 0) {
        real = pow(a.data[0], b);
        if (!isnan(real)) {
            return quaternion_from_re(real);
        }
    }
    value = complex_lpows(complex_make(a.data[0], radius), b);
    return quaternion_complex_slice(a, value, radius);
}

static inline quaternion
quaternion_rpow(quaternion a, quaternion b)
{
    if (vec4_is_real(b)) {
        return quaternion_lpows(a, b.data[0]);
    }
    return quaternion_exp(quaternion_mul(quaternion_log(a), b));
}

static inline quaternion
quaternion_pow(quaternion a, quaternion b)
{
    return quaternion_rpow(a, b);
}

static inline quaternion
quaternion_lpow(quaternion a, quaternion b)
{
    return quaternion_exp(quaternion_mul(b, quaternion_log(a)));
}

static inline quaternion
quaternion_radd(double a, quaternion b) {
  return vec4_add(quaternion_from_re(a), b);
}

static inline quaternion
quaternion_ladd(quaternion a, double b) {
  return vec4_add(a, quaternion_from_re(b));
}

static inline quaternion
quaternion_rsub(double a, quaternion b) {
  return vec4_sub(quaternion_from_re(a), b);
}

static inline quaternion
quaternion_lsub(quaternion a, double b) {
  return vec4_sub(a, quaternion_from_re(b));
}

static inline quaternion
quaternion_rmul(double a, quaternion b) {
  return quaternion_mul(quaternion_from_re(a), b);
}

static inline quaternion
quaternion_lmul(quaternion a, double b) {
  return quaternion_mul(a, quaternion_from_re(b));
}

static inline quaternion
quaternion_rdiv(double a, quaternion b) {
  return quaternion_div(quaternion_from_re(a), b);
}

static inline quaternion
quaternion_ldiv(quaternion a, double b) {
  return quaternion_div(a, quaternion_from_re(b));
}

static inline quaternion
quaternion_rpows(double a, quaternion b) {
  return quaternion_rpow(quaternion_from_re(a), b);
}

static inline quaternion
quaternion_mix(quaternion a, quaternion b, double t)
{
    quaternion delta = quaternion_mul(quaternion_inv(a), b);
    return quaternion_mul(a, quaternion_lpows(delta, t));
}

static inline void
quaternion_minpoly(quaternion a, double result[3])
{
    result[0] = quaternion_norm2(a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/*
 * Eigenvalues x0 -/+ sqrt(-|v|^2) of the quaternion, as complex numbers
 * stored in a vec2 (data[0] = real part, data[1] = imaginary part).
 * Always a conjugate pair; real quaternions give a double eigenvalue.
 */
static inline void
quaternion_eivals(quaternion a, complex result[2])
{
    double vn = quaternion_vector_norm(a);
    result[0] = complex_make(a.data[0], -vn);
    result[1] = complex_make(a.data[0], +vn);
}

/* Common algebra API. */
static inline quaternion
quaternion_zero(void)
{
    return vec4_zero();
}

static inline int
quaternion_is_real(quaternion a)
{
    return vec4_is_real(a);
}

static inline quaternion
quaternion_pos(quaternion a)
{
    return vec4_pos(a);
}

static inline quaternion
quaternion_neg(quaternion a)
{
    return vec4_neg(a);
}

static inline quaternion
quaternion_add(quaternion a, quaternion b)
{
    return vec4_add(a, b);
}

static inline quaternion
quaternion_sub(quaternion a, quaternion b)
{
    return vec4_sub(a, b);
}

static inline int
quaternion_equal(quaternion a, quaternion b)
{
    return vec4_all_eq(a, b);
}

static inline int
quaternion_not_equal(quaternion a, quaternion b)
{
    return !quaternion_equal(a, b);
}

static inline quaternion
quaternion_sqrt(quaternion a)
{
    return quaternion_apply(a, complex_sqrt);
}

static inline quaternion
quaternion_cos(quaternion a)
{
    return quaternion_apply(a, complex_cos);
}

static inline quaternion
quaternion_sin(quaternion a)
{
    return quaternion_apply(a, complex_sin);
}

static inline quaternion
quaternion_tan(quaternion a)
{
    return quaternion_apply(a, complex_tan);
}

static inline quaternion
quaternion_acos(quaternion a)
{
    return quaternion_apply(a, complex_acos);
}

static inline quaternion
quaternion_asin(quaternion a)
{
    return quaternion_apply(a, complex_asin);
}

static inline quaternion
quaternion_atan(quaternion a)
{
    return quaternion_apply(a, complex_atan);
}

static inline quaternion
quaternion_cosh(quaternion a)
{
    return quaternion_apply(a, complex_cosh);
}

static inline quaternion
quaternion_sinh(quaternion a)
{
    return quaternion_apply(a, complex_sinh);
}

static inline quaternion
quaternion_tanh(quaternion a)
{
    return quaternion_apply(a, complex_tanh);
}

static inline quaternion
quaternion_acosh(quaternion a)
{
    return quaternion_apply(a, complex_acosh);
}

static inline quaternion
quaternion_asinh(quaternion a)
{
    return quaternion_apply(a, complex_asinh);
}

static inline quaternion
quaternion_atanh(quaternion a)
{
    return quaternion_apply(a, complex_atanh);
}



static inline quaternion
quaternion_expm1(quaternion a)
{
    return quaternion_apply(a, complex_expm1);
}

static inline quaternion
quaternion_log1p(quaternion a)
{
    return quaternion_apply(a, complex_log1p);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
quaternion_lmat(quaternion a, double result[ALGPY_QUATERNION_N*ALGPY_QUATERNION_N])
{
    quaternion basis = {{0}};
    quaternion product;
    int row, col;
    for (col = 0; col < ALGPY_QUATERNION_N; col++) {
        basis.data[col] = 1;
        product = quaternion_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_QUATERNION_N; row++) {
            result[row * ALGPY_QUATERNION_N + col] = product.data[row];
        }
    }
}

static inline void
quaternion_rmat(quaternion a, double result[ALGPY_QUATERNION_N*ALGPY_QUATERNION_N])
{
    quaternion basis = {{0}};
    quaternion product;
    int row, col;
    for (row = 0; row < ALGPY_QUATERNION_N; row++) {
        basis.data[row] = 1;
        product = quaternion_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_QUATERNION_N; col++) {
            result[row * ALGPY_QUATERNION_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_QUATERNION_H */
