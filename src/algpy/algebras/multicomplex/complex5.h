#ifndef ALGPY_COMPLEX5_H
#define ALGPY_COMPLEX5_H

#include <string.h>
#include <stdint.h>
#include "../../vectorspaces/vec32.h"
#include "complex4.h"

#define ALGPY_COMPLEX5_N ALGPY_VEC32_N

typedef vec32 complex5;

static inline complex5
complex5_make(const double coeffs[ALGPY_COMPLEX5_N])
{
    complex5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

/* Write a = u + v*i4 and evaluate f(u +/- i0*v) in complex4. Reconstruct
 * u' = (p+q)/2 and v' = (p-q)/(2*i0), as in complex4_apply.
 * Principal branches are selected independently in each complex component. */
static inline complex5
complex5_apply(complex5 a, complex4 (*f)(complex4))
{
    complex4 p = complex4_make((const double[]){
        a.data[0] - a.data[17],
        a.data[1] + a.data[16],
        a.data[2] - a.data[19],
        a.data[3] + a.data[18],
        a.data[4] - a.data[21],
        a.data[5] + a.data[20],
        a.data[6] - a.data[23],
        a.data[7] + a.data[22],
        a.data[8] - a.data[25],
        a.data[9] + a.data[24],
        a.data[10] - a.data[27],
        a.data[11] + a.data[26],
        a.data[12] - a.data[29],
        a.data[13] + a.data[28],
        a.data[14] - a.data[31],
        a.data[15] + a.data[30]});
    complex4 q = complex4_make((const double[]){
        a.data[0] + a.data[17],
        a.data[1] - a.data[16],
        a.data[2] + a.data[19],
        a.data[3] - a.data[18],
        a.data[4] + a.data[21],
        a.data[5] - a.data[20],
        a.data[6] + a.data[23],
        a.data[7] - a.data[22],
        a.data[8] + a.data[25],
        a.data[9] - a.data[24],
        a.data[10] + a.data[27],
        a.data[11] - a.data[26],
        a.data[12] + a.data[29],
        a.data[13] - a.data[28],
        a.data[14] + a.data[31],
        a.data[15] - a.data[30]});
    p = f(p);
    q = f(q);
    return complex5_make((const double[]){
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[2] + q.data[2]) / 2,
        (p.data[3] + q.data[3]) / 2,
        (p.data[4] + q.data[4]) / 2,
        (p.data[5] + q.data[5]) / 2,
        (p.data[6] + q.data[6]) / 2,
        (p.data[7] + q.data[7]) / 2,
        (p.data[8] + q.data[8]) / 2,
        (p.data[9] + q.data[9]) / 2,
        (p.data[10] + q.data[10]) / 2,
        (p.data[11] + q.data[11]) / 2,
        (p.data[12] + q.data[12]) / 2,
        (p.data[13] + q.data[13]) / 2,
        (p.data[14] + q.data[14]) / 2,
        (p.data[15] + q.data[15]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (q.data[0] - p.data[0]) / 2,
        (p.data[3] - q.data[3]) / 2,
        (q.data[2] - p.data[2]) / 2,
        (p.data[5] - q.data[5]) / 2,
        (q.data[4] - p.data[4]) / 2,
        (p.data[7] - q.data[7]) / 2,
        (q.data[6] - p.data[6]) / 2,
        (p.data[9] - q.data[9]) / 2,
        (q.data[8] - p.data[8]) / 2,
        (p.data[11] - q.data[11]) / 2,
        (q.data[10] - p.data[10]) / 2,
        (p.data[13] - q.data[13]) / 2,
        (q.data[12] - p.data[12]) / 2,
        (p.data[15] - q.data[15]) / 2,
        (q.data[14] - p.data[14]) / 2});
}

static inline complex5
complex5_from_re(double a)
{
    complex5 result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
complex5_to_re(complex5 a)
{
    return a.data[0];
}

static inline complex5
complex5_zero(void)
{
    return vec32_zero();
}

static inline complex5
complex5_one(void)
{
    return complex5_from_re(1);
}

static inline int
complex5_is_real(complex5 a)
{
    return vec32_is_real(a);
}

static inline complex5
complex5_pos(complex5 a)
{
    return vec32_pos(a);
}

static inline complex5
complex5_neg(complex5 a)
{
    return vec32_neg(a);
}

static inline complex5
complex5_re(complex5 a)
{
    return complex5_from_re(a.data[0]);
}

static inline complex5
complex5_im(complex5 a)
{
    return vec32_sub(a, complex5_re(a));
}

static inline complex5
complex5_conjugate(complex5 a)
{
    return complex5_make((const double[]){a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15], -a.data[16], -a.data[17], -a.data[18], -a.data[19], -a.data[20], -a.data[21], -a.data[22], -a.data[23], -a.data[24], -a.data[25], -a.data[26], -a.data[27], -a.data[28], -a.data[29], -a.data[30], -a.data[31]});
}

static inline complex5
complex5_add(complex5 a, complex5 b)
{
    return vec32_add(a, b);
}

static inline complex5
complex5_sub(complex5 a, complex5 b)
{
    return vec32_sub(a, b);
}

static inline int
complex5_equal(complex5 a, complex5 b)
{
    return vec32_all_eq(a, b);
}

static inline int
complex5_not_equal(complex5 a, complex5 b)
{
    return !complex5_equal(a, b);
}

static inline complex5
complex5_mul(complex5 a, complex5 b)
{
    if (complex5_is_real(b)) {
        return vec32_scale(a, b.data[0]);
    }
    if (complex5_is_real(a)) {
        return vec32_scale(b, a.data[0]);
    }
    complex4 u = complex4_make(a.data);
    complex4 v = complex4_make(a.data + 16);
    complex4 s = complex4_make(b.data);
    complex4 t = complex4_make(b.data + 16);
    complex4 lo = complex4_sub(complex4_mul(u, s), complex4_mul(v, t));
    complex4 hi = complex4_add(complex4_mul(u, t), complex4_mul(v, s));
    complex5 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_COMPLEX5_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

/* Recursive norm: (u + v*i)*(u - v*i) = u*u + v*v. */
static inline complex4
complex5_norm2(complex5 a)
{
    complex4 u = complex4_make(a.data);
    complex4 v = complex4_make(a.data + 16);
    return complex4_add(complex4_mul(u, u), complex4_mul(v, v));
}

static inline complex3
complex5_norm4(complex5 a)
{
    return complex4_norm2(complex5_norm2(a));
}

static inline complex2
complex5_norm8(complex5 a)
{
    return complex3_norm2(complex5_norm4(a));
}

static inline complex1
complex5_norm16(complex5 a)
{
    return complex2_norm2(complex5_norm8(a));
}

static inline double
complex5_norm32(complex5 a)
{
    return complex1_norm2(complex5_norm16(a));
}

static inline double
complex5_norm(complex5 a)
{
    return sqrt(sqrt(sqrt(sqrt(sqrt(complex5_norm32(a))))));
}

static inline complex5
complex5_inv(complex5 a)
{
    complex4 ninv = complex4_inv(complex5_norm2(a));
    complex4 u = complex4_make(a.data);
    complex4 v = complex4_make(a.data + 16);
    complex4 lo = complex4_mul(u, ninv);
    complex4 hi = complex4_mul(complex4_neg(v), ninv);
    complex5 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_COMPLEX5_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

static inline complex5
complex5_div(complex5 a, complex5 b)
{
    if (complex5_is_real(b)) {
        return vec32_div_scalar(a, b.data[0]);
    }
    return complex5_mul(a, complex5_inv(b));
}

static inline complex5
complex5_square(complex5 a)
{
    return complex5_mul(a, a);
}

static inline complex5
complex5_exp(complex5 a)
{
    return complex5_apply(a, complex4_exp);
}

static inline complex5
complex5_log(complex5 a)
{
    return complex5_apply(a, complex4_log);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline complex5
complex5_ipow(complex5 a, int64_t b)
{
    complex5 result = complex5_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return complex5_mul(a, a);
    }
    if (b == -1) {
        return complex5_inv(a);
    }
    if (b < 0) {
        a = complex5_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = complex5_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = complex5_mul(a, a);
        }
    }
    return result;
}

static inline complex5
complex5_rpow(complex5 a, complex5 b)
{
    return complex5_exp(complex5_mul(complex5_log(a), b));
}

static inline complex5
complex5_pow(complex5 a, complex5 b)
{
    return complex5_rpow(a, b);
}

static inline complex5
complex5_lpow(complex5 a, complex5 b)
{
    return complex5_rpow(a, b);
}

static inline complex5
complex5_sqrt(complex5 a)
{
    return complex5_apply(a, complex4_sqrt);
}

static inline complex5
complex5_cos(complex5 a)
{
    return complex5_apply(a, complex4_cos);
}

static inline complex5
complex5_sin(complex5 a)
{
    return complex5_apply(a, complex4_sin);
}

static inline complex5
complex5_tan(complex5 a)
{
    return complex5_apply(a, complex4_tan);
}

static inline complex5
complex5_acos(complex5 a)
{
    return complex5_apply(a, complex4_acos);
}

static inline complex5
complex5_asin(complex5 a)
{
    return complex5_apply(a, complex4_asin);
}

static inline complex5
complex5_atan(complex5 a)
{
    return complex5_apply(a, complex4_atan);
}

static inline complex5
complex5_cosh(complex5 a)
{
    return complex5_apply(a, complex4_cosh);
}

static inline complex5
complex5_sinh(complex5 a)
{
    return complex5_apply(a, complex4_sinh);
}

static inline complex5
complex5_tanh(complex5 a)
{
    return complex5_apply(a, complex4_tanh);
}

static inline complex5
complex5_acosh(complex5 a)
{
    return complex5_apply(a, complex4_acosh);
}

static inline complex5
complex5_asinh(complex5 a)
{
    return complex5_apply(a, complex4_asinh);
}

static inline complex5
complex5_atanh(complex5 a)
{
    return complex5_apply(a, complex4_atanh);
}

static inline complex5
complex5_expm1(complex5 a)
{
    return complex5_apply(a, complex4_expm1);
}

static inline complex5
complex5_log1p(complex5 a)
{
    return complex5_apply(a, complex4_log1p);
}

static inline complex5
complex5_radd(double a, complex5 b)
{
    return complex5_add(complex5_from_re(a), b);
}

static inline complex5
complex5_ladd(complex5 a, double b)
{
    return complex5_add(a, complex5_from_re(b));
}

static inline complex5
complex5_rsub(double a, complex5 b)
{
    return complex5_sub(complex5_from_re(a), b);
}

static inline complex5
complex5_lsub(complex5 a, double b)
{
    return complex5_sub(a, complex5_from_re(b));
}

static inline complex5
complex5_rmul(double a, complex5 b)
{
    return vec32_scale(b, a);
}

static inline complex5
complex5_lmul(complex5 a, double b)
{
    return vec32_scale(a, b);
}

static inline complex5
complex5_rdiv(double a, complex5 b)
{
    return complex5_div(complex5_from_re(a), b);
}

static inline complex5
complex5_ldiv(complex5 a, double b)
{
    return vec32_div_scalar(a, b);
}

static inline complex5
complex5_rpows(double a, complex5 b)
{
    return complex5_rpow(complex5_from_re(a), b);
}

static inline complex5
complex5_lpows(complex5 a, double b)
{
    return complex5_rpow(a, complex5_from_re(b));
}

static inline complex5
complex5_mix(complex5 a, complex5 b, double t)
{
    complex5 delta = complex5_mul(complex5_inv(a), b);
    return complex5_mul(a, complex5_lpows(delta, t));
}

static inline void
complex5_basis(complex5 result[32])
{
    int i, j;
    for (i = 0; i < ALGPY_COMPLEX5_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX5_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
complex5_gens(complex5 result[5])
{
    result[0] = complex5_make((const double[ALGPY_COMPLEX5_N]){[1] = 1});
    result[1] = complex5_make((const double[ALGPY_COMPLEX5_N]){[2] = 1});
    result[2] = complex5_make((const double[ALGPY_COMPLEX5_N]){[4] = 1});
    result[3] = complex5_make((const double[ALGPY_COMPLEX5_N]){[8] = 1});
    result[4] = complex5_make((const double[ALGPY_COMPLEX5_N]){[16] = 1});
}

static inline void
complex5_coeffs(complex5 result[1024])
{
    complex5 basis[ALGPY_COMPLEX5_N];
    int i, j;
    complex5_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX5_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX5_N; j++) {
            result[i*ALGPY_COMPLEX5_N+j] = complex5_mul(basis[i], basis[j]);
        }
    }
}

static inline void
complex5_minpoly(complex5 a, double result[33])
{
    complex4 p = complex4_make((const double[]){
        a.data[0] - a.data[17],
        a.data[1] + a.data[16],
        a.data[2] - a.data[19],
        a.data[3] + a.data[18],
        a.data[4] - a.data[21],
        a.data[5] + a.data[20],
        a.data[6] - a.data[23],
        a.data[7] + a.data[22],
        a.data[8] - a.data[25],
        a.data[9] + a.data[24],
        a.data[10] - a.data[27],
        a.data[11] + a.data[26],
        a.data[12] - a.data[29],
        a.data[13] + a.data[28],
        a.data[14] - a.data[31],
        a.data[15] + a.data[30]});
    complex4 q = complex4_make((const double[]){
        a.data[0] + a.data[17],
        a.data[1] - a.data[16],
        a.data[2] + a.data[19],
        a.data[3] - a.data[18],
        a.data[4] + a.data[21],
        a.data[5] - a.data[20],
        a.data[6] + a.data[23],
        a.data[7] - a.data[22],
        a.data[8] + a.data[25],
        a.data[9] - a.data[24],
        a.data[10] + a.data[27],
        a.data[11] - a.data[26],
        a.data[12] + a.data[29],
        a.data[13] - a.data[28],
        a.data[14] + a.data[31],
        a.data[15] - a.data[30]});
    double left[17], right[17];
    int i, j;
    complex4_minpoly(p, left);
    complex4_minpoly(q, right);
    for (i = 0; i <= 32; i++) result[i] = 0;
    for (i = 0; i <= 16; i++) {
        for (j = 0; j <= 16; j++) result[i+j] += left[i]*right[j];
    }
}

static inline void
complex5_eivals(complex5 a, complex result[32])
{
    complex4 p = complex4_make((const double[]){
        a.data[0] - a.data[17],
        a.data[1] + a.data[16],
        a.data[2] - a.data[19],
        a.data[3] + a.data[18],
        a.data[4] - a.data[21],
        a.data[5] + a.data[20],
        a.data[6] - a.data[23],
        a.data[7] + a.data[22],
        a.data[8] - a.data[25],
        a.data[9] + a.data[24],
        a.data[10] - a.data[27],
        a.data[11] + a.data[26],
        a.data[12] - a.data[29],
        a.data[13] + a.data[28],
        a.data[14] - a.data[31],
        a.data[15] + a.data[30]});
    complex4 q = complex4_make((const double[]){
        a.data[0] + a.data[17],
        a.data[1] - a.data[16],
        a.data[2] + a.data[19],
        a.data[3] - a.data[18],
        a.data[4] + a.data[21],
        a.data[5] - a.data[20],
        a.data[6] + a.data[23],
        a.data[7] - a.data[22],
        a.data[8] + a.data[25],
        a.data[9] - a.data[24],
        a.data[10] + a.data[27],
        a.data[11] - a.data[26],
        a.data[12] + a.data[29],
        a.data[13] - a.data[28],
        a.data[14] + a.data[31],
        a.data[15] - a.data[30]});
    complex left[16], right[16];
    complex4_eivals(p, left);
    complex4_eivals(q, right);
    /* Preserve the sign-bit ordering of the generators. */
    result[0] = left[0];
    result[1] = right[1];
    result[2] = left[2];
    result[3] = right[3];
    result[4] = left[4];
    result[5] = right[5];
    result[6] = left[6];
    result[7] = right[7];
    result[8] = left[8];
    result[9] = right[9];
    result[10] = left[10];
    result[11] = right[11];
    result[12] = left[12];
    result[13] = right[13];
    result[14] = left[14];
    result[15] = right[15];
    result[16] = right[0];
    result[17] = left[1];
    result[18] = right[2];
    result[19] = left[3];
    result[20] = right[4];
    result[21] = left[5];
    result[22] = right[6];
    result[23] = left[7];
    result[24] = right[8];
    result[25] = left[9];
    result[26] = right[10];
    result[27] = left[11];
    result[28] = right[12];
    result[29] = left[13];
    result[30] = right[14];
    result[31] = left[15];
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
complex5_lmat(complex5 a, double result[ALGPY_COMPLEX5_N*ALGPY_COMPLEX5_N])
{
    complex5 basis = {{0}};
    complex5 product;
    int row, col;
    for (col = 0; col < ALGPY_COMPLEX5_N; col++) {
        basis.data[col] = 1;
        product = complex5_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_COMPLEX5_N; row++) {
            result[row * ALGPY_COMPLEX5_N + col] = product.data[row];
        }
    }
}

static inline void
complex5_rmat(complex5 a, double result[ALGPY_COMPLEX5_N*ALGPY_COMPLEX5_N])
{
    complex5 basis = {{0}};
    complex5 product;
    int row, col;
    for (row = 0; row < ALGPY_COMPLEX5_N; row++) {
        basis.data[row] = 1;
        product = complex5_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_COMPLEX5_N; col++) {
            result[row * ALGPY_COMPLEX5_N + col] = product.data[col];
        }
    }
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline complex5
complex5_grade0(complex5 a)
{
    const double *x = a.data;
    return complex5_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline complex5
complex5_grade1(complex5 a)
{
    const double *x = a.data;
    return complex5_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline complex5
complex5_grade2(complex5 a)
{
    const double *x = a.data;
    return complex5_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline complex5
complex5_grade3(complex5 a)
{
    const double *x = a.data;
    return complex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline complex5
complex5_grade4(complex5 a)
{
    const double *x = a.data;
    return complex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline complex5
complex5_grade5(complex5 a)
{
    const double *x = a.data;
    return complex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

#endif /* ALGPY_COMPLEX5_H */
