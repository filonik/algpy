/*
 * Core octonion arithmetic with fixed generator squares -1.
 *
 * Pure C with no dependency on Python or NumPy, so that it can be reused and
 * tested independently of the wrapper code in `_octonion.c`.
 */
#ifndef ALGPY_OCTONION_H
#define ALGPY_OCTONION_H

#include <stdint.h>
#include "../../fields/complex.h"
#include "quaternion.h"
#include "../../vectorspaces/vec2.h"
#include "../../vectorspaces/vec8.h"

#define ALGPY_OCTONION_N ALGPY_VEC8_N

typedef vec8 octonion;

static inline octonion
octonion_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
    octonion result = {{x0, x1, x2, x3, x4, x5, x6, x7}};
    return result;
}

static inline octonion
octonion_re(octonion a)
{
    return octonion_make(a.data[0], 0, 0, 0, 0, 0, 0, 0);
}

static inline octonion
octonion_im(octonion a)
{
    return octonion_make(0, a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
}

static inline octonion
octonion_from_re(double a)
{
    return octonion_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
octonion_to_re(octonion a)
{
    return a.data[0];
}

static inline octonion
octonion_one(void)
{
    return octonion_make(1, 0, 0, 0, 0, 0, 0, 0);
}

/* Fixed-parameter Cayley-Dickson product of quaternion pairs. */
static inline octonion
octonion_mul(octonion a, octonion b)
{
    if (vec8_is_real(b)) {
        return vec8_scale(a, b.data[0]);
    }
    if (vec8_is_real(a)) {
        return vec8_scale(b, a.data[0]);
    }
    quaternion a0 = {{a.data[0], a.data[1], a.data[2], a.data[3]}};
    quaternion a1 = {{a.data[4], a.data[5], a.data[6], a.data[7]}};
    quaternion b0 = {{b.data[0], b.data[1], b.data[2], b.data[3]}};
    quaternion b1 = {{b.data[4], b.data[5], b.data[6], b.data[7]}};
    quaternion left = vec4_sub(quaternion_mul(a0, b0),
        quaternion_mul(quaternion_conjugate(b1), a1));
    quaternion right = vec4_add(quaternion_mul(b1, a0),
        quaternion_mul(a1, quaternion_conjugate(b0)));
    return octonion_make(
        left.data[0], left.data[1], left.data[2], left.data[3],
        right.data[0], right.data[1], right.data[2], right.data[3]);
}

static inline void
octonion_basis(octonion result[ALGPY_OCTONION_N])
{
    int i, j;
    for (i = 0; i < ALGPY_OCTONION_N; i++) {
            for (j = 0; j < ALGPY_OCTONION_N; j++) {
                result[i].data[j] = i == j ? 1 : 0;
            }
    }
}

static inline void
octonion_coeffs(octonion result[ALGPY_OCTONION_N * ALGPY_OCTONION_N])
{
    octonion basis[ALGPY_OCTONION_N];
    int i, j;
    octonion_basis(basis);
    for (i = 0; i < ALGPY_OCTONION_N; i++) {
            for (j = 0; j < ALGPY_OCTONION_N; j++) {
                result[i * ALGPY_OCTONION_N + j] =
                    octonion_mul(basis[i], basis[j]);
            }
    }
}

static inline double
octonion_norm2(octonion a)
{
    return a.data[0]*a.data[0]
        + a.data[1]*a.data[1]
        + a.data[2]*a.data[2]
        + a.data[3]*a.data[3]
        + a.data[4]*a.data[4]
        + a.data[5]*a.data[5]
        + a.data[6]*a.data[6]
        + a.data[7]*a.data[7];
}

static inline double
octonion_norm(octonion a)
{
    return sqrt(octonion_norm2(a));
}

static inline double
octonion_vector_norm2(octonion a)
{
    return a.data[1]*a.data[1]
        + a.data[2]*a.data[2]
        + a.data[3]*a.data[3]
        + a.data[4]*a.data[4]
        + a.data[5]*a.data[5]
        + a.data[6]*a.data[6]
        + a.data[7]*a.data[7];
}

static inline double
octonion_vector_norm(octonion a)
{
    double r0 = hypot(a.data[1], a.data[2]);
    double r1 = hypot(a.data[3], a.data[4]);
    double r2 = hypot(a.data[5], a.data[6]);
    return hypot(hypot(r0, r1), hypot(r2, a.data[7]));
}

static inline octonion
octonion_conjugate(octonion a)
{
    return octonion_make(
        a.data[0], -a.data[1], -a.data[2], -a.data[3],
        -a.data[4], -a.data[5], -a.data[6], -a.data[7]);
}

static inline octonion
octonion_inv(octonion a)
{
    return vec8_div_scalar(octonion_conjugate(a), octonion_norm2(a));
}

static inline octonion
octonion_div(octonion a, octonion b)
{
    if (vec8_is_real(b)) {
        return vec8_div_scalar(a, b.data[0]);
    }
    return octonion_mul(a, octonion_inv(b));
}

/* Lift a complex value along the input's imaginary direction. Real inputs
 * use the first imaginary basis vector to select the principal branch. */
static inline octonion
octonion_complex_slice(octonion a, complex value, double radius)
{
    octonion result = {{0}};
    int i;
    result.data[0] = value.data[0];
    if (radius == 0) {
        result.data[1] = value.data[1];
        return result;
    }
    for (i = 1; i < ALGPY_OCTONION_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
    }
    return result;
}

static inline octonion
octonion_apply(octonion a, complex (*f)(complex))
{
    double radius = octonion_vector_norm(a);
    complex value = f(complex_make(a.data[0], radius));
    return octonion_complex_slice(a, value, radius);
}

static inline octonion
octonion_exp(octonion a)
{
    return octonion_apply(a, complex_exp);
}

/* Principal logarithm; for negative reals the pure part points along i. */
static inline octonion
octonion_log(octonion a)
{
    return octonion_apply(a, complex_log);
}

static inline octonion
octonion_square(octonion a)
{
    return octonion_make(
        a.data[0]*a.data[0] - octonion_vector_norm2(a), 2*a.data[0]*a.data[1], 2*a.data[0]*a.data[2], 2*a.data[0]*a.data[3],
        2*a.data[0]*a.data[4], 2*a.data[0]*a.data[5], 2*a.data[0]*a.data[6], 2*a.data[0]*a.data[7]);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline octonion
octonion_ipow(octonion a, int64_t b)
{
    octonion result = octonion_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return octonion_mul(a, a);
    }
    if (b == -1) {
        return octonion_inv(a);
    }
    if (b < 0) {
        a = octonion_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = octonion_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = octonion_mul(a, a);
        }
    }
    return result;
}

static inline octonion
octonion_lpows(octonion a, double b)
{
    double radius, real;
    complex value;
    if (b == 0) {
        return octonion_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return octonion_square(a);
    }
    radius = octonion_vector_norm(a);
    if (radius == 0) {
        real = pow(a.data[0], b);
        if (!isnan(real)) {
            return octonion_from_re(real);
        }
    }
    value = complex_lpows(complex_make(a.data[0], radius), b);
    return octonion_complex_slice(a, value, radius);
}

static inline octonion
octonion_rpow(octonion a, octonion b)
{
    if (vec8_is_real(b)) {
        return octonion_lpows(a, b.data[0]);
    }
    return octonion_exp(octonion_mul(octonion_log(a), b));
}

static inline octonion
octonion_pow(octonion a, octonion b)
{
    return octonion_rpow(a, b);
}

static inline octonion
octonion_lpow(octonion a, octonion b)
{
    return octonion_exp(octonion_mul(b, octonion_log(a)));
}

static inline octonion
octonion_radd(double a, octonion b) {
  return vec8_add(octonion_from_re(a), b);
}

static inline octonion
octonion_ladd(octonion a, double b) {
  return vec8_add(a, octonion_from_re(b));
}

static inline octonion
octonion_rsub(double a, octonion b) {
  return vec8_sub(octonion_from_re(a), b);
}

static inline octonion
octonion_lsub(octonion a, double b) {
  return vec8_sub(a, octonion_from_re(b));
}

static inline octonion
octonion_rmul(double a, octonion b) {
  return octonion_mul(octonion_from_re(a), b);
}

static inline octonion
octonion_lmul(octonion a, double b) {
  return octonion_mul(a, octonion_from_re(b));
}

static inline octonion
octonion_rdiv(double a, octonion b) {
  return octonion_div(octonion_from_re(a), b);
}

static inline octonion
octonion_ldiv(octonion a, double b) {
  return octonion_div(a, octonion_from_re(b));
}

static inline octonion
octonion_rpows(double a, octonion b) {
  return octonion_rpow(octonion_from_re(a), b);
}

static inline octonion
octonion_mix(octonion a, octonion b, double t)
{
    octonion delta = octonion_mul(octonion_inv(a), b);
    return octonion_mul(a, octonion_lpows(delta, t));
}

static inline void
octonion_minpoly(octonion a, double result[3])
{
    result[0] = octonion_norm2(a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/*
 * Eigenvalues x0 -/+ sqrt(-|v|^2) of the octonion, as complex numbers
 * stored in a vec2 (data[0] = real part, data[1] = imaginary part).
 * Always a conjugate pair; real octonions give a double eigenvalue.
 */
static inline void
octonion_eivals(octonion a, complex result[2])
{
    double vn = octonion_vector_norm(a);
    result[0] = complex_make(a.data[0], -vn);
    result[1] = complex_make(a.data[0], +vn);
}

/* Common algebra API. */
static inline octonion
octonion_zero(void)
{
    return vec8_zero();
}

static inline int
octonion_is_real(octonion a)
{
    return vec8_is_real(a);
}

static inline octonion
octonion_pos(octonion a)
{
    return vec8_pos(a);
}

static inline octonion
octonion_neg(octonion a)
{
    return vec8_neg(a);
}

static inline octonion
octonion_add(octonion a, octonion b)
{
    return vec8_add(a, b);
}

static inline octonion
octonion_sub(octonion a, octonion b)
{
    return vec8_sub(a, b);
}

static inline int
octonion_equal(octonion a, octonion b)
{
    return vec8_all_eq(a, b);
}

static inline int
octonion_not_equal(octonion a, octonion b)
{
    return !octonion_equal(a, b);
}

static inline octonion
octonion_sqrt(octonion a)
{
    return octonion_apply(a, complex_sqrt);
}

static inline octonion
octonion_cos(octonion a)
{
    return octonion_apply(a, complex_cos);
}

static inline octonion
octonion_sin(octonion a)
{
    return octonion_apply(a, complex_sin);
}

static inline octonion
octonion_tan(octonion a)
{
    return octonion_apply(a, complex_tan);
}

static inline octonion
octonion_acos(octonion a)
{
    return octonion_apply(a, complex_acos);
}

static inline octonion
octonion_asin(octonion a)
{
    return octonion_apply(a, complex_asin);
}

static inline octonion
octonion_atan(octonion a)
{
    return octonion_apply(a, complex_atan);
}

static inline octonion
octonion_cosh(octonion a)
{
    return octonion_apply(a, complex_cosh);
}

static inline octonion
octonion_sinh(octonion a)
{
    return octonion_apply(a, complex_sinh);
}

static inline octonion
octonion_tanh(octonion a)
{
    return octonion_apply(a, complex_tanh);
}

static inline octonion
octonion_acosh(octonion a)
{
    return octonion_apply(a, complex_acosh);
}

static inline octonion
octonion_asinh(octonion a)
{
    return octonion_apply(a, complex_asinh);
}

static inline octonion
octonion_atanh(octonion a)
{
    return octonion_apply(a, complex_atanh);
}



static inline octonion
octonion_expm1(octonion a)
{
    return octonion_apply(a, complex_expm1);
}

static inline octonion
octonion_log1p(octonion a)
{
    return octonion_apply(a, complex_log1p);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
octonion_lmat(octonion a, double result[ALGPY_OCTONION_N*ALGPY_OCTONION_N])
{
    octonion basis = {{0}};
    octonion product;
    int row, col;
    for (col = 0; col < ALGPY_OCTONION_N; col++) {
        basis.data[col] = 1;
        product = octonion_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_OCTONION_N; row++) {
            result[row * ALGPY_OCTONION_N + col] = product.data[row];
        }
    }
}

static inline void
octonion_rmat(octonion a, double result[ALGPY_OCTONION_N*ALGPY_OCTONION_N])
{
    octonion basis = {{0}};
    octonion product;
    int row, col;
    for (row = 0; row < ALGPY_OCTONION_N; row++) {
        basis.data[row] = 1;
        product = octonion_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_OCTONION_N; col++) {
            result[row * ALGPY_OCTONION_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_OCTONION_H */
