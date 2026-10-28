/*
 * Core quaternion arithmetic.
 *
 * Pure C with no dependency on Python or NumPy, so that it can be reused and
 * tested independently of the wrapper code in `_quaternion.c`.
 */
#ifndef ALGPY_QUATERNION_H
#define ALGPY_QUATERNION_H

#include "../../fields/complex.h"
#include "../../vectorspaces/vec2.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_QUATERNION_N ALGPY_VEC4_N

typedef vec4 quaternion;

static inline quaternion
quaternion_make(double w, double x, double y, double z)
{
    quaternion q = {{w, x, y, z}};
    return q;
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
    return a.data[0] * a.data[0] + a.data[1] * a.data[1] + a.data[2] * a.data[2] + a.data[3] * a.data[3];
}

static inline double
quaternion_norm(quaternion a)
{
    return sqrt(quaternion_norm2(a));
}

static inline double
quaternion_vector_norm2(quaternion a)
{
    return a.data[1] * a.data[1] + a.data[2] * a.data[2] + a.data[3] * a.data[3];
}

static inline double
quaternion_vector_norm(quaternion a)
{
    return sqrt(quaternion_vector_norm2(a));
}

static inline quaternion
quaternion_conjugate(quaternion a)
{
    return quaternion_make(a.data[0], -a.data[1], -a.data[2], -a.data[3]);
}

static inline quaternion
quaternion_inv(quaternion a)
{
    double n = quaternion_norm2(a);
    return quaternion_make(a.data[0] / n, -a.data[1] / n, -a.data[2] / n, -a.data[3] / n);
}

static inline quaternion
quaternion_div(quaternion a, quaternion b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return quaternion_mul(a, quaternion_inv(b));
}

static inline quaternion
quaternion_exp(quaternion a)
{
    double vn = quaternion_vector_norm(a);
    double ew = exp(a.data[0]);
    double s;
    if (vn == 0) {
        return quaternion_make(ew, 0, 0, 0);
    }
    s = ew * sin(vn) / vn;
    return quaternion_make(ew * cos(vn), s * a.data[1], s * a.data[2],
                           s * a.data[3]);
}

/* Principal logarithm; for negative reals the pure part points along i. */
static inline quaternion
quaternion_log(quaternion a)
{
    double vn = quaternion_vector_norm(a);
    double f;
    if (vn == 0) {
        if (a.data[0] < 0) {
            return quaternion_make(log(-a.data[0]), M_PI, 0, 0);
        }
        return quaternion_make(log(a.data[0]), 0, 0, 0);
    }
    f = atan2(vn, a.data[0]) / vn;
    return quaternion_make(log(sqrt(quaternion_norm2(a))),
                           f * a.data[1], f * a.data[2], f * a.data[3]);
}

static inline quaternion
quaternion_pow(quaternion a, quaternion b)
{
    return quaternion_exp(quaternion_mul(quaternion_log(a), b));
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
quaternion_rpow(double a, quaternion b) {
  return quaternion_pow(quaternion_from_re(a), b);
}

static inline quaternion
quaternion_lpow(quaternion a, double b) {
  return quaternion_pow(a, quaternion_from_re(b));
}

static inline quaternion
quaternion_mix(quaternion a, quaternion b, double t)
{
    quaternion delta = quaternion_mul(quaternion_inv(a), b);
    return quaternion_mul(a, quaternion_lpow(delta, t));
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

/* Common algebra API. Unfinished operations return NaN components. */
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
    /* TODO: implement quaternion_sqrt. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_cos(quaternion a)
{
    /* TODO: implement quaternion_cos. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_sin(quaternion a)
{
    /* TODO: implement quaternion_sin. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_tan(quaternion a)
{
    /* TODO: implement quaternion_tan. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_acos(quaternion a)
{
    /* TODO: implement quaternion_acos. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_asin(quaternion a)
{
    /* TODO: implement quaternion_asin. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_atan(quaternion a)
{
    /* TODO: implement quaternion_atan. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_cosh(quaternion a)
{
    /* TODO: implement quaternion_cosh. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_sinh(quaternion a)
{
    /* TODO: implement quaternion_sinh. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_tanh(quaternion a)
{
    /* TODO: implement quaternion_tanh. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_acosh(quaternion a)
{
    /* TODO: implement quaternion_acosh. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_asinh(quaternion a)
{
    /* TODO: implement quaternion_asinh. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_atanh(quaternion a)
{
    /* TODO: implement quaternion_atanh. */
    return quaternion_make(NAN, NAN, NAN, NAN);
}

static inline quaternion
quaternion_square(quaternion a)
{
    return quaternion_mul(a, a);
}

static inline quaternion
quaternion_expm1(quaternion a)
{
    return quaternion_lsub(quaternion_exp(a), 1);
}

static inline quaternion
quaternion_log1p(quaternion a)
{
    return quaternion_log(quaternion_ladd(a, 1));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
quaternion_lmat(quaternion a, double result[ALGPY_QUATERNION_N*ALGPY_QUATERNION_N])
{
    quaternion coeffs[ALGPY_QUATERNION_N * ALGPY_QUATERNION_N];
    int row, col, i;
    quaternion_coeffs(coeffs);
    for (row = 0; row < ALGPY_QUATERNION_N; row++) {
        for (col = 0; col < ALGPY_QUATERNION_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_QUATERNION_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_QUATERNION_N + col].data[row];
            }
            result[row * ALGPY_QUATERNION_N + col] = value;
        }
    }
}

static inline void
quaternion_rmat(quaternion a, double result[ALGPY_QUATERNION_N*ALGPY_QUATERNION_N])
{
    quaternion coeffs[ALGPY_QUATERNION_N * ALGPY_QUATERNION_N];
    int row, col, i;
    quaternion_coeffs(coeffs);
    for (row = 0; row < ALGPY_QUATERNION_N; row++) {
        for (col = 0; col < ALGPY_QUATERNION_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_QUATERNION_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_QUATERNION_N + i].data[col];
            }
            result[row * ALGPY_QUATERNION_N + col] = value;
        }
    }
}


#endif /* ALGPY_QUATERNION_H */
