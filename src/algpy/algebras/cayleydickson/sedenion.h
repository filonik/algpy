/*
 * Core sedenion arithmetic with fixed generator squares -1.
 *
 * Pure C with no dependency on Python or NumPy, so that it can be reused and
 * tested independently of the wrapper code in `_sedenion.c`.
 */
#ifndef ALGPY_SEDENION_H
#define ALGPY_SEDENION_H

#include <string.h>
#include <stdint.h>
#include "../../fields/complex.h"
#include "octonion.h"
#include "../../vectorspaces/vec2.h"
#include "../../vectorspaces/vec16.h"

#define ALGPY_SEDENION_N ALGPY_VEC16_N

typedef vec16 sedenion;

static inline sedenion
sedenion_make(const double coeffs[ALGPY_SEDENION_N])
{
    sedenion result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline sedenion
sedenion_re(sedenion a)
{
    return sedenion_make((const double[]){a.data[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline sedenion
sedenion_im(sedenion a)
{
    return sedenion_make((const double[]){
        0, a.data[1], a.data[2], a.data[3],
        a.data[4], a.data[5], a.data[6], a.data[7],
        a.data[8], a.data[9], a.data[10], a.data[11],
        a.data[12], a.data[13], a.data[14], a.data[15]});
}

static inline sedenion
sedenion_from_re(double a)
{
    return sedenion_make((const double[ALGPY_SEDENION_N]){a});
}

static inline double
sedenion_to_re(sedenion a)
{
    return a.data[0];
}

static inline sedenion
sedenion_one(void)
{
    return sedenion_make((const double[ALGPY_SEDENION_N]){[0] = 1});
}

/* Fixed-parameter Cayley-Dickson product of octonion pairs. */
static inline sedenion
sedenion_mul(sedenion a, sedenion b)
{
    if (vec16_is_real(b)) {
        return vec16_scale(a, b.data[0]);
    }
    if (vec16_is_real(a)) {
        return vec16_scale(b, a.data[0]);
    }
    octonion a0 = {{a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]}};
    octonion a1 = {{a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]}};
    octonion b0 = {{b.data[0], b.data[1], b.data[2], b.data[3], b.data[4], b.data[5], b.data[6], b.data[7]}};
    octonion b1 = {{b.data[8], b.data[9], b.data[10], b.data[11], b.data[12], b.data[13], b.data[14], b.data[15]}};
    octonion left = vec8_sub(octonion_mul(a0, b0),
        octonion_mul(octonion_conjugate(b1), a1));
    octonion right = vec8_add(octonion_mul(b1, a0),
        octonion_mul(a1, octonion_conjugate(b0)));
    sedenion result;
    memcpy(result.data, left.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_SEDENION_N/2, right.data, sizeof(result.data)/2);
    return result;
}

static inline void
sedenion_basis(sedenion result[ALGPY_SEDENION_N])
{
    int i, j;
    for (i = 0; i < ALGPY_SEDENION_N; i++) {
            for (j = 0; j < ALGPY_SEDENION_N; j++) {
                result[i].data[j] = i == j ? 1 : 0;
            }
    }
}

static inline void
sedenion_gens(sedenion result[4])
{
    result[0] = sedenion_make((const double[ALGPY_SEDENION_N]){[1] = 1});
    result[1] = sedenion_make((const double[ALGPY_SEDENION_N]){[2] = 1});
    result[2] = sedenion_make((const double[ALGPY_SEDENION_N]){[4] = 1});
    result[3] = sedenion_make((const double[ALGPY_SEDENION_N]){[8] = 1});
}

static inline void
sedenion_coeffs(sedenion result[ALGPY_SEDENION_N * ALGPY_SEDENION_N])
{
    sedenion basis[ALGPY_SEDENION_N];
    int i, j;
    sedenion_basis(basis);
    for (i = 0; i < ALGPY_SEDENION_N; i++) {
            for (j = 0; j < ALGPY_SEDENION_N; j++) {
                result[i * ALGPY_SEDENION_N + j] =
                    sedenion_mul(basis[i], basis[j]);
            }
    }
}

static inline double
sedenion_norm2(sedenion a)
{
    return a.data[0]*a.data[0]
        + a.data[1]*a.data[1]
        + a.data[2]*a.data[2]
        + a.data[3]*a.data[3]
        + a.data[4]*a.data[4]
        + a.data[5]*a.data[5]
        + a.data[6]*a.data[6]
        + a.data[7]*a.data[7]
        + a.data[8]*a.data[8]
        + a.data[9]*a.data[9]
        + a.data[10]*a.data[10]
        + a.data[11]*a.data[11]
        + a.data[12]*a.data[12]
        + a.data[13]*a.data[13]
        + a.data[14]*a.data[14]
        + a.data[15]*a.data[15];
}

static inline double
sedenion_norm(sedenion a)
{
    return sqrt(sedenion_norm2(a));
}

static inline double
sedenion_vector_norm2(sedenion a)
{
    return a.data[1]*a.data[1]
        + a.data[2]*a.data[2]
        + a.data[3]*a.data[3]
        + a.data[4]*a.data[4]
        + a.data[5]*a.data[5]
        + a.data[6]*a.data[6]
        + a.data[7]*a.data[7]
        + a.data[8]*a.data[8]
        + a.data[9]*a.data[9]
        + a.data[10]*a.data[10]
        + a.data[11]*a.data[11]
        + a.data[12]*a.data[12]
        + a.data[13]*a.data[13]
        + a.data[14]*a.data[14]
        + a.data[15]*a.data[15];
}

static inline double
sedenion_vector_norm(sedenion a)
{
    double r0 = hypot(a.data[1], a.data[2]);
    double r1 = hypot(a.data[3], a.data[4]);
    double r2 = hypot(a.data[5], a.data[6]);
    double r3 = hypot(a.data[7], a.data[8]);
    double r4 = hypot(a.data[9], a.data[10]);
    double r5 = hypot(a.data[11], a.data[12]);
    double r6 = hypot(a.data[13], a.data[14]);
    return hypot(hypot(hypot(r0, r1), hypot(r2, r3)), hypot(hypot(r4, r5), hypot(r6, a.data[15])));
}

static inline sedenion
sedenion_conjugate(sedenion a)
{
    return sedenion_make((const double[]){
        a.data[0], -a.data[1], -a.data[2], -a.data[3],
        -a.data[4], -a.data[5], -a.data[6], -a.data[7],
        -a.data[8], -a.data[9], -a.data[10], -a.data[11],
        -a.data[12], -a.data[13], -a.data[14], -a.data[15]});
}

static inline sedenion
sedenion_inv(sedenion a)
{
    return vec16_div_scalar(sedenion_conjugate(a), sedenion_norm2(a));
}

static inline sedenion
sedenion_div(sedenion a, sedenion b)
{
    if (vec16_is_real(b)) {
        return vec16_div_scalar(a, b.data[0]);
    }
    return sedenion_mul(a, sedenion_inv(b));
}

/* Lift a complex value along the input's imaginary direction. Real inputs
 * use the first imaginary basis vector to select the principal branch. */
static inline sedenion
sedenion_complex_slice(sedenion a, complex value, double radius)
{
    sedenion result = {{0}};
    int i;
    result.data[0] = value.data[0];
    if (radius == 0) {
        result.data[1] = value.data[1];
        return result;
    }
    for (i = 1; i < ALGPY_SEDENION_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
    }
    return result;
}

static inline sedenion
sedenion_apply(sedenion a, complex (*f)(complex))
{
    double radius = sedenion_vector_norm(a);
    complex value = f(complex_make((const double[]){a.data[0], radius}));
    return sedenion_complex_slice(a, value, radius);
}

static inline sedenion
sedenion_exp(sedenion a)
{
    return sedenion_apply(a, complex_exp);
}

/* Principal logarithm; for negative reals the pure part points along i. */
static inline sedenion
sedenion_log(sedenion a)
{
    return sedenion_apply(a, complex_log);
}

static inline sedenion
sedenion_square(sedenion a)
{
    return sedenion_make((const double[]){
        a.data[0]*a.data[0] - sedenion_vector_norm2(a), 2*a.data[0]*a.data[1], 2*a.data[0]*a.data[2], 2*a.data[0]*a.data[3],
        2*a.data[0]*a.data[4], 2*a.data[0]*a.data[5], 2*a.data[0]*a.data[6], 2*a.data[0]*a.data[7],
        2*a.data[0]*a.data[8], 2*a.data[0]*a.data[9], 2*a.data[0]*a.data[10], 2*a.data[0]*a.data[11],
        2*a.data[0]*a.data[12], 2*a.data[0]*a.data[13], 2*a.data[0]*a.data[14], 2*a.data[0]*a.data[15]});
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline sedenion
sedenion_ipow(sedenion a, int64_t b)
{
    sedenion result = sedenion_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return sedenion_mul(a, a);
    }
    if (b == -1) {
        return sedenion_inv(a);
    }
    if (b < 0) {
        a = sedenion_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = sedenion_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = sedenion_mul(a, a);
        }
    }
    return result;
}

static inline sedenion
sedenion_lpows(sedenion a, double b)
{
    double radius, real;
    complex value;
    if (b == 0) {
        return sedenion_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return sedenion_square(a);
    }
    radius = sedenion_vector_norm(a);
    if (radius == 0) {
        real = pow(a.data[0], b);
        if (!isnan(real)) {
            return sedenion_from_re(real);
        }
    }
    value = complex_lpows(complex_make((const double[]){a.data[0], radius}), b);
    return sedenion_complex_slice(a, value, radius);
}

static inline sedenion
sedenion_rpow(sedenion a, sedenion b)
{
    if (vec16_is_real(b)) {
        return sedenion_lpows(a, b.data[0]);
    }
    return sedenion_exp(sedenion_mul(sedenion_log(a), b));
}

static inline sedenion
sedenion_pow(sedenion a, sedenion b)
{
    return sedenion_rpow(a, b);
}

static inline sedenion
sedenion_lpow(sedenion a, sedenion b)
{
    return sedenion_exp(sedenion_mul(b, sedenion_log(a)));
}

static inline sedenion
sedenion_radd(double a, sedenion b) {
  return vec16_add(sedenion_from_re(a), b);
}

static inline sedenion
sedenion_ladd(sedenion a, double b) {
  return vec16_add(a, sedenion_from_re(b));
}

static inline sedenion
sedenion_rsub(double a, sedenion b) {
  return vec16_sub(sedenion_from_re(a), b);
}

static inline sedenion
sedenion_lsub(sedenion a, double b) {
  return vec16_sub(a, sedenion_from_re(b));
}

static inline sedenion
sedenion_rmul(double a, sedenion b) {
  return sedenion_mul(sedenion_from_re(a), b);
}

static inline sedenion
sedenion_lmul(sedenion a, double b) {
  return sedenion_mul(a, sedenion_from_re(b));
}

static inline sedenion
sedenion_rdiv(double a, sedenion b) {
  return sedenion_div(sedenion_from_re(a), b);
}

static inline sedenion
sedenion_ldiv(sedenion a, double b) {
  return sedenion_div(a, sedenion_from_re(b));
}

static inline sedenion
sedenion_rpows(double a, sedenion b) {
  return sedenion_rpow(sedenion_from_re(a), b);
}

static inline sedenion
sedenion_mix(sedenion a, sedenion b, double t)
{
    sedenion delta = sedenion_mul(sedenion_inv(a), b);
    return sedenion_mul(a, sedenion_lpows(delta, t));
}

static inline void
sedenion_minpoly(sedenion a, double result[3])
{
    result[0] = sedenion_norm2(a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/*
 * Eigenvalues x0 -/+ sqrt(-|v|^2) of the sedenion, as complex numbers
 * stored in a vec2 (data[0] = real part, data[1] = imaginary part).
 * Always a conjugate pair; real sedenions give a double eigenvalue.
 */
static inline void
sedenion_eivals(sedenion a, complex result[2])
{
    double vn = sedenion_vector_norm(a);
    result[0] = complex_make((const double[]){a.data[0], -vn});
    result[1] = complex_make((const double[]){a.data[0], +vn});
}

/* Common algebra API. */
static inline sedenion
sedenion_zero(void)
{
    return vec16_zero();
}

static inline int
sedenion_is_real(sedenion a)
{
    return vec16_is_real(a);
}

static inline sedenion
sedenion_pos(sedenion a)
{
    return vec16_pos(a);
}

static inline sedenion
sedenion_neg(sedenion a)
{
    return vec16_neg(a);
}

static inline sedenion
sedenion_add(sedenion a, sedenion b)
{
    return vec16_add(a, b);
}

static inline sedenion
sedenion_sub(sedenion a, sedenion b)
{
    return vec16_sub(a, b);
}

static inline int
sedenion_equal(sedenion a, sedenion b)
{
    return vec16_all_eq(a, b);
}

static inline int
sedenion_not_equal(sedenion a, sedenion b)
{
    return !sedenion_equal(a, b);
}

static inline sedenion
sedenion_sqrt(sedenion a)
{
    return sedenion_apply(a, complex_sqrt);
}

static inline sedenion
sedenion_cos(sedenion a)
{
    return sedenion_apply(a, complex_cos);
}

static inline sedenion
sedenion_sin(sedenion a)
{
    return sedenion_apply(a, complex_sin);
}

static inline sedenion
sedenion_tan(sedenion a)
{
    return sedenion_apply(a, complex_tan);
}

static inline sedenion
sedenion_acos(sedenion a)
{
    return sedenion_apply(a, complex_acos);
}

static inline sedenion
sedenion_asin(sedenion a)
{
    return sedenion_apply(a, complex_asin);
}

static inline sedenion
sedenion_atan(sedenion a)
{
    return sedenion_apply(a, complex_atan);
}

static inline sedenion
sedenion_cosh(sedenion a)
{
    return sedenion_apply(a, complex_cosh);
}

static inline sedenion
sedenion_sinh(sedenion a)
{
    return sedenion_apply(a, complex_sinh);
}

static inline sedenion
sedenion_tanh(sedenion a)
{
    return sedenion_apply(a, complex_tanh);
}

static inline sedenion
sedenion_acosh(sedenion a)
{
    return sedenion_apply(a, complex_acosh);
}

static inline sedenion
sedenion_asinh(sedenion a)
{
    return sedenion_apply(a, complex_asinh);
}

static inline sedenion
sedenion_atanh(sedenion a)
{
    return sedenion_apply(a, complex_atanh);
}



static inline sedenion
sedenion_expm1(sedenion a)
{
    return sedenion_apply(a, complex_expm1);
}

static inline sedenion
sedenion_log1p(sedenion a)
{
    return sedenion_apply(a, complex_log1p);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
sedenion_lmat(sedenion a, double result[ALGPY_SEDENION_N*ALGPY_SEDENION_N])
{
    sedenion basis = {{0}};
    sedenion product;
    int row, col;
    for (col = 0; col < ALGPY_SEDENION_N; col++) {
        basis.data[col] = 1;
        product = sedenion_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_SEDENION_N; row++) {
            result[row * ALGPY_SEDENION_N + col] = product.data[row];
        }
    }
}

static inline void
sedenion_rmat(sedenion a, double result[ALGPY_SEDENION_N*ALGPY_SEDENION_N])
{
    sedenion basis = {{0}};
    sedenion product;
    int row, col;
    for (row = 0; row < ALGPY_SEDENION_N; row++) {
        basis.data[row] = 1;
        product = sedenion_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_SEDENION_N; col++) {
            result[row * ALGPY_SEDENION_N + col] = product.data[col];
        }
    }
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline sedenion
sedenion_grade0(sedenion a)
{
    const double *x = a.data;
    return sedenion_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline sedenion
sedenion_grade1(sedenion a)
{
    const double *x = a.data;
    return sedenion_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0});
}

static inline sedenion
sedenion_grade2(sedenion a)
{
    const double *x = a.data;
    return sedenion_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0});
}

static inline sedenion
sedenion_grade3(sedenion a)
{
    const double *x = a.data;
    return sedenion_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0});
}

static inline sedenion
sedenion_grade4(sedenion a)
{
    const double *x = a.data;
    return sedenion_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15]});
}

#endif /* ALGPY_SEDENION_H */
