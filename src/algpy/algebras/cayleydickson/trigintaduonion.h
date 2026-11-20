/*
 * Core trigintaduonion arithmetic with fixed generator squares -1.
 *
 * Pure C with no dependency on Python or NumPy, so that it can be reused and
 * tested independently of the wrapper code in `_trigintaduonion.c`.
 */
#ifndef ALGPY_TRIGINTADUONION_H
#define ALGPY_TRIGINTADUONION_H

#include <string.h>
#include <stdint.h>
#include "../../fields/complex.h"
#include "sedenion.h"
#include "../../vectorspaces/vec2.h"
#include "../../vectorspaces/vec32.h"

#define ALGPY_TRIGINTADUONION_N ALGPY_VEC32_N

typedef vec32 trigintaduonion;

static inline trigintaduonion
trigintaduonion_make(const double coeffs[ALGPY_TRIGINTADUONION_N])
{
    trigintaduonion result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline trigintaduonion
trigintaduonion_re(trigintaduonion a)
{
    trigintaduonion result = {{0}};
    result.data[0] = a.data[0];
    return result;
}

static inline trigintaduonion
trigintaduonion_im(trigintaduonion a)
{
    return trigintaduonion_make((const double[]){0, a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15], a.data[16], a.data[17], a.data[18], a.data[19], a.data[20], a.data[21], a.data[22], a.data[23], a.data[24], a.data[25], a.data[26], a.data[27], a.data[28], a.data[29], a.data[30], a.data[31]});
}

static inline trigintaduonion
trigintaduonion_from_re(double a)
{
    trigintaduonion result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
trigintaduonion_to_re(trigintaduonion a)
{
    return a.data[0];
}

static inline trigintaduonion
trigintaduonion_one(void)
{
    trigintaduonion result = {{0}};
    result.data[0] = 1;
    return result;
}

/* Fixed-parameter Cayley-Dickson product of sedenion pairs. */
static inline trigintaduonion
trigintaduonion_mul(trigintaduonion a, trigintaduonion b)
{
    if (vec32_is_real(b)) {
        return vec32_scale(a, b.data[0]);
    }
    if (vec32_is_real(a)) {
        return vec32_scale(b, a.data[0]);
    }
    sedenion a0 = {{a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]}};
    sedenion a1 = {{a.data[16], a.data[17], a.data[18], a.data[19], a.data[20], a.data[21], a.data[22], a.data[23], a.data[24], a.data[25], a.data[26], a.data[27], a.data[28], a.data[29], a.data[30], a.data[31]}};
    sedenion b0 = {{b.data[0], b.data[1], b.data[2], b.data[3], b.data[4], b.data[5], b.data[6], b.data[7], b.data[8], b.data[9], b.data[10], b.data[11], b.data[12], b.data[13], b.data[14], b.data[15]}};
    sedenion b1 = {{b.data[16], b.data[17], b.data[18], b.data[19], b.data[20], b.data[21], b.data[22], b.data[23], b.data[24], b.data[25], b.data[26], b.data[27], b.data[28], b.data[29], b.data[30], b.data[31]}};
    sedenion left = vec16_sub(sedenion_mul(a0, b0),
        sedenion_mul(sedenion_conjugate(b1), a1));
    sedenion right = vec16_add(sedenion_mul(b1, a0),
        sedenion_mul(a1, sedenion_conjugate(b0)));
    trigintaduonion result;
    memcpy(result.data, left.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_TRIGINTADUONION_N/2, right.data, sizeof(result.data)/2);
    return result;
}

static inline void
trigintaduonion_basis(trigintaduonion result[ALGPY_TRIGINTADUONION_N])
{
    int i, j;
    for (i = 0; i < ALGPY_TRIGINTADUONION_N; i++) {
        for (j = 0; j < ALGPY_TRIGINTADUONION_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
trigintaduonion_gens(trigintaduonion result[5])
{
    result[0] = trigintaduonion_make((const double[ALGPY_TRIGINTADUONION_N]){[1] = 1});
    result[1] = trigintaduonion_make((const double[ALGPY_TRIGINTADUONION_N]){[2] = 1});
    result[2] = trigintaduonion_make((const double[ALGPY_TRIGINTADUONION_N]){[4] = 1});
    result[3] = trigintaduonion_make((const double[ALGPY_TRIGINTADUONION_N]){[8] = 1});
    result[4] = trigintaduonion_make((const double[ALGPY_TRIGINTADUONION_N]){[16] = 1});
}

static inline void
trigintaduonion_coeffs(trigintaduonion result[ALGPY_TRIGINTADUONION_N * ALGPY_TRIGINTADUONION_N])
{
    trigintaduonion basis[ALGPY_TRIGINTADUONION_N];
    int i, j;
    trigintaduonion_basis(basis);
    for (i = 0; i < ALGPY_TRIGINTADUONION_N; i++) {
        for (j = 0; j < ALGPY_TRIGINTADUONION_N; j++) {
            result[i * ALGPY_TRIGINTADUONION_N + j] =
                trigintaduonion_mul(basis[i], basis[j]);
        }
    }
}

static inline double
trigintaduonion_norm2(trigintaduonion a)
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
        + a.data[15]*a.data[15]
        + a.data[16]*a.data[16]
        + a.data[17]*a.data[17]
        + a.data[18]*a.data[18]
        + a.data[19]*a.data[19]
        + a.data[20]*a.data[20]
        + a.data[21]*a.data[21]
        + a.data[22]*a.data[22]
        + a.data[23]*a.data[23]
        + a.data[24]*a.data[24]
        + a.data[25]*a.data[25]
        + a.data[26]*a.data[26]
        + a.data[27]*a.data[27]
        + a.data[28]*a.data[28]
        + a.data[29]*a.data[29]
        + a.data[30]*a.data[30]
        + a.data[31]*a.data[31];
}

static inline double
trigintaduonion_norm(trigintaduonion a)
{
    return sqrt(trigintaduonion_norm2(a));
}

static inline double
trigintaduonion_vector_norm2(trigintaduonion a)
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
        + a.data[15]*a.data[15]
        + a.data[16]*a.data[16]
        + a.data[17]*a.data[17]
        + a.data[18]*a.data[18]
        + a.data[19]*a.data[19]
        + a.data[20]*a.data[20]
        + a.data[21]*a.data[21]
        + a.data[22]*a.data[22]
        + a.data[23]*a.data[23]
        + a.data[24]*a.data[24]
        + a.data[25]*a.data[25]
        + a.data[26]*a.data[26]
        + a.data[27]*a.data[27]
        + a.data[28]*a.data[28]
        + a.data[29]*a.data[29]
        + a.data[30]*a.data[30]
        + a.data[31]*a.data[31];
}

static inline double
trigintaduonion_vector_norm(trigintaduonion a)
{
    double r0 = hypot(a.data[1], a.data[2]);
    double r1 = hypot(a.data[3], a.data[4]);
    double r2 = hypot(a.data[5], a.data[6]);
    double r3 = hypot(a.data[7], a.data[8]);
    double r4 = hypot(a.data[9], a.data[10]);
    double r5 = hypot(a.data[11], a.data[12]);
    double r6 = hypot(a.data[13], a.data[14]);
    double r7 = hypot(a.data[15], a.data[16]);
    double r8 = hypot(a.data[17], a.data[18]);
    double r9 = hypot(a.data[19], a.data[20]);
    double r10 = hypot(a.data[21], a.data[22]);
    double r11 = hypot(a.data[23], a.data[24]);
    double r12 = hypot(a.data[25], a.data[26]);
    double r13 = hypot(a.data[27], a.data[28]);
    double r14 = hypot(a.data[29], a.data[30]);
    return hypot(
        hypot(hypot(hypot(r0, r1), hypot(r2, r3)),
              hypot(hypot(r4, r5), hypot(r6, r7))),
        hypot(hypot(hypot(r8, r9), hypot(r10, r11)),
              hypot(hypot(r12, r13), hypot(r14, a.data[31]))));
}

static inline trigintaduonion
trigintaduonion_conjugate(trigintaduonion a)
{
    return trigintaduonion_make((const double[]){a.data[0], -a.data[1], -a.data[2], -a.data[3], -a.data[4], -a.data[5], -a.data[6], -a.data[7], -a.data[8], -a.data[9], -a.data[10], -a.data[11], -a.data[12], -a.data[13], -a.data[14], -a.data[15], -a.data[16], -a.data[17], -a.data[18], -a.data[19], -a.data[20], -a.data[21], -a.data[22], -a.data[23], -a.data[24], -a.data[25], -a.data[26], -a.data[27], -a.data[28], -a.data[29], -a.data[30], -a.data[31]});
}

static inline trigintaduonion
trigintaduonion_inv(trigintaduonion a)
{
    return vec32_div_scalar(trigintaduonion_conjugate(a), trigintaduonion_norm2(a));
}

static inline trigintaduonion
trigintaduonion_div(trigintaduonion a, trigintaduonion b)
{
    if (vec32_is_real(b)) {
        return vec32_div_scalar(a, b.data[0]);
    }
    return trigintaduonion_mul(a, trigintaduonion_inv(b));
}

/* Lift a complex value along the input's imaginary direction. Real inputs
 * use the first imaginary basis vector to select the principal branch. */
static inline trigintaduonion
trigintaduonion_complex_slice(trigintaduonion a, complex value, double radius)
{
    trigintaduonion result = {{0}};
    int i;
    result.data[0] = value.data[0];
    if (radius == 0) {
        result.data[1] = value.data[1];
        return result;
    }
    for (i = 1; i < ALGPY_TRIGINTADUONION_N; i++) {
        result.data[i] = a.data[i] == 0 ? 0 : (a.data[i]/radius)*value.data[1];
    }
    return result;
}

static inline trigintaduonion
trigintaduonion_apply(trigintaduonion a, complex (*f)(complex))
{
    double radius = trigintaduonion_vector_norm(a);
    complex value = f(complex_make((const double[]){a.data[0], radius}));
    return trigintaduonion_complex_slice(a, value, radius);
}

static inline trigintaduonion
trigintaduonion_exp(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_exp);
}

/* Principal logarithm; for negative reals the pure part points along i. */
static inline trigintaduonion
trigintaduonion_log(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_log);
}

static inline trigintaduonion
trigintaduonion_square(trigintaduonion a)
{
    return trigintaduonion_make((const double[]){a.data[0]*a.data[0] - trigintaduonion_vector_norm2(a), 2*a.data[0]*a.data[1], 2*a.data[0]*a.data[2], 2*a.data[0]*a.data[3], 2*a.data[0]*a.data[4], 2*a.data[0]*a.data[5], 2*a.data[0]*a.data[6], 2*a.data[0]*a.data[7], 2*a.data[0]*a.data[8], 2*a.data[0]*a.data[9], 2*a.data[0]*a.data[10], 2*a.data[0]*a.data[11], 2*a.data[0]*a.data[12], 2*a.data[0]*a.data[13], 2*a.data[0]*a.data[14], 2*a.data[0]*a.data[15], 2*a.data[0]*a.data[16], 2*a.data[0]*a.data[17], 2*a.data[0]*a.data[18], 2*a.data[0]*a.data[19], 2*a.data[0]*a.data[20], 2*a.data[0]*a.data[21], 2*a.data[0]*a.data[22], 2*a.data[0]*a.data[23], 2*a.data[0]*a.data[24], 2*a.data[0]*a.data[25], 2*a.data[0]*a.data[26], 2*a.data[0]*a.data[27], 2*a.data[0]*a.data[28], 2*a.data[0]*a.data[29], 2*a.data[0]*a.data[30], 2*a.data[0]*a.data[31]});
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline trigintaduonion
trigintaduonion_ipow(trigintaduonion a, int64_t b)
{
    trigintaduonion result = trigintaduonion_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return trigintaduonion_mul(a, a);
    }
    if (b == -1) {
        return trigintaduonion_inv(a);
    }
    if (b < 0) {
        a = trigintaduonion_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = trigintaduonion_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = trigintaduonion_mul(a, a);
        }
    }
    return result;
}

static inline trigintaduonion
trigintaduonion_lpows(trigintaduonion a, double b)
{
    double radius, real;
    complex value;
    radius = trigintaduonion_vector_norm(a);
    if (radius == 0) {
        real = pow(a.data[0], b);
        if (!isnan(real)) {
            return trigintaduonion_from_re(real);
        }
    }
    value = complex_lpows(complex_make((const double[]){a.data[0], radius}), b);
    return trigintaduonion_complex_slice(a, value, radius);
}

static inline trigintaduonion
trigintaduonion_rpow(trigintaduonion a, trigintaduonion b)
{
    if (vec32_is_real(b)) {
        return trigintaduonion_lpows(a, b.data[0]);
    }
    return trigintaduonion_exp(trigintaduonion_mul(trigintaduonion_log(a), b));
}

static inline trigintaduonion
trigintaduonion_pow(trigintaduonion a, trigintaduonion b)
{
    return trigintaduonion_rpow(a, b);
}

static inline trigintaduonion
trigintaduonion_lpow(trigintaduonion a, trigintaduonion b)
{
    return trigintaduonion_exp(trigintaduonion_mul(b, trigintaduonion_log(a)));
}

static inline trigintaduonion
trigintaduonion_radd(double a, trigintaduonion b)
{
    return vec32_add(trigintaduonion_from_re(a), b);
}

static inline trigintaduonion
trigintaduonion_ladd(trigintaduonion a, double b)
{
    return vec32_add(a, trigintaduonion_from_re(b));
}

static inline trigintaduonion
trigintaduonion_rsub(double a, trigintaduonion b)
{
    return vec32_sub(trigintaduonion_from_re(a), b);
}

static inline trigintaduonion
trigintaduonion_lsub(trigintaduonion a, double b)
{
    return vec32_sub(a, trigintaduonion_from_re(b));
}

static inline trigintaduonion
trigintaduonion_rmul(double a, trigintaduonion b)
{
    return trigintaduonion_mul(trigintaduonion_from_re(a), b);
}

static inline trigintaduonion
trigintaduonion_lmul(trigintaduonion a, double b)
{
    return trigintaduonion_mul(a, trigintaduonion_from_re(b));
}

static inline trigintaduonion
trigintaduonion_rdiv(double a, trigintaduonion b)
{
    return trigintaduonion_div(trigintaduonion_from_re(a), b);
}

static inline trigintaduonion
trigintaduonion_ldiv(trigintaduonion a, double b)
{
    return trigintaduonion_div(a, trigintaduonion_from_re(b));
}

static inline trigintaduonion
trigintaduonion_rpows(double a, trigintaduonion b)
{
    return trigintaduonion_rpow(trigintaduonion_from_re(a), b);
}

static inline trigintaduonion
trigintaduonion_mix(trigintaduonion a, trigintaduonion b, double t)
{
    trigintaduonion delta = trigintaduonion_mul(trigintaduonion_inv(a), b);
    return trigintaduonion_mul(a, trigintaduonion_lpows(delta, t));
}

static inline void
trigintaduonion_minpoly(trigintaduonion a, double result[3])
{
    result[0] = trigintaduonion_norm2(a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/*
 * Eigenvalues x0 -/+ sqrt(-|v|^2) of the trigintaduonion, as complex numbers
 * stored in a vec2 (data[0] = real part, data[1] = imaginary part).
 * Always a conjugate pair; real trigintaduonions give a double eigenvalue.
 */
static inline void
trigintaduonion_eivals(trigintaduonion a, complex result[2])
{
    double vn = trigintaduonion_vector_norm(a);
    result[0] = complex_make((const double[]){a.data[0], -vn});
    result[1] = complex_make((const double[]){a.data[0], +vn});
}

/* Common algebra API. */
static inline trigintaduonion
trigintaduonion_zero(void)
{
    return vec32_zero();
}

static inline int
trigintaduonion_is_real(trigintaduonion a)
{
    return vec32_is_real(a);
}

static inline trigintaduonion
trigintaduonion_pos(trigintaduonion a)
{
    return vec32_pos(a);
}

static inline trigintaduonion
trigintaduonion_neg(trigintaduonion a)
{
    return vec32_neg(a);
}

static inline trigintaduonion
trigintaduonion_add(trigintaduonion a, trigintaduonion b)
{
    return vec32_add(a, b);
}

static inline trigintaduonion
trigintaduonion_sub(trigintaduonion a, trigintaduonion b)
{
    return vec32_sub(a, b);
}

static inline int
trigintaduonion_equal(trigintaduonion a, trigintaduonion b)
{
    return vec32_all_eq(a, b);
}

static inline int
trigintaduonion_not_equal(trigintaduonion a, trigintaduonion b)
{
    return !trigintaduonion_equal(a, b);
}

static inline trigintaduonion
trigintaduonion_sqrt(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_sqrt);
}

static inline trigintaduonion
trigintaduonion_cos(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_cos);
}

static inline trigintaduonion
trigintaduonion_sin(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_sin);
}

static inline trigintaduonion
trigintaduonion_tan(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_tan);
}

static inline trigintaduonion
trigintaduonion_acos(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_acos);
}

static inline trigintaduonion
trigintaduonion_asin(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_asin);
}

static inline trigintaduonion
trigintaduonion_atan(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_atan);
}

static inline trigintaduonion
trigintaduonion_cosh(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_cosh);
}

static inline trigintaduonion
trigintaduonion_sinh(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_sinh);
}

static inline trigintaduonion
trigintaduonion_tanh(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_tanh);
}

static inline trigintaduonion
trigintaduonion_acosh(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_acosh);
}

static inline trigintaduonion
trigintaduonion_asinh(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_asinh);
}

static inline trigintaduonion
trigintaduonion_atanh(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_atanh);
}

static inline trigintaduonion
trigintaduonion_expm1(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_expm1);
}

static inline trigintaduonion
trigintaduonion_log1p(trigintaduonion a)
{
    return trigintaduonion_apply(a, complex_log1p);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
trigintaduonion_lmat(trigintaduonion a, double result[ALGPY_TRIGINTADUONION_N*ALGPY_TRIGINTADUONION_N])
{
    trigintaduonion basis = {{0}};
    trigintaduonion product;
    int row, col;
    for (col = 0; col < ALGPY_TRIGINTADUONION_N; col++) {
        basis.data[col] = 1;
        product = trigintaduonion_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_TRIGINTADUONION_N; row++) {
            result[row * ALGPY_TRIGINTADUONION_N + col] = product.data[row];
        }
    }
}

static inline void
trigintaduonion_rmat(trigintaduonion a, double result[ALGPY_TRIGINTADUONION_N*ALGPY_TRIGINTADUONION_N])
{
    trigintaduonion basis = {{0}};
    trigintaduonion product;
    int row, col;
    for (row = 0; row < ALGPY_TRIGINTADUONION_N; row++) {
        basis.data[row] = 1;
        product = trigintaduonion_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_TRIGINTADUONION_N; col++) {
            result[row * ALGPY_TRIGINTADUONION_N + col] = product.data[col];
        }
    }
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline trigintaduonion
trigintaduonion_grade0(trigintaduonion a)
{
    const double *x = a.data;
    return trigintaduonion_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline trigintaduonion
trigintaduonion_grade1(trigintaduonion a)
{
    const double *x = a.data;
    return trigintaduonion_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline trigintaduonion
trigintaduonion_grade2(trigintaduonion a)
{
    const double *x = a.data;
    return trigintaduonion_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline trigintaduonion
trigintaduonion_grade3(trigintaduonion a)
{
    const double *x = a.data;
    return trigintaduonion_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline trigintaduonion
trigintaduonion_grade4(trigintaduonion a)
{
    const double *x = a.data;
    return trigintaduonion_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline trigintaduonion
trigintaduonion_grade5(trigintaduonion a)
{
    const double *x = a.data;
    return trigintaduonion_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

#endif /* ALGPY_TRIGINTADUONION_H */
