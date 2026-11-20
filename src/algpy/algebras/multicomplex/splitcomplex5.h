#ifndef ALGPY_SPLITCOMPLEX5_H
#define ALGPY_SPLITCOMPLEX5_H

#include <string.h>
#include <stdint.h>
#include "../../vectorspaces/vec32.h"
#include "splitcomplex4.h"

#define ALGPY_SPLITCOMPLEX5_N ALGPY_VEC32_N

typedef vec32 splitcomplex5;

static inline splitcomplex5
splitcomplex5_make(const double coeffs[ALGPY_SPLITCOMPLEX5_N])
{
    splitcomplex5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

/* For a = u + v*e with e*e = +1, apply f to u+v and u-v. */
static inline splitcomplex5
splitcomplex5_apply(splitcomplex5 a, splitcomplex4 (*f)(splitcomplex4))
{
    splitcomplex4 p = splitcomplex4_make((const double[]){
        a.data[0] + a.data[16],
        a.data[1] + a.data[17],
        a.data[2] + a.data[18],
        a.data[3] + a.data[19],
        a.data[4] + a.data[20],
        a.data[5] + a.data[21],
        a.data[6] + a.data[22],
        a.data[7] + a.data[23],
        a.data[8] + a.data[24],
        a.data[9] + a.data[25],
        a.data[10] + a.data[26],
        a.data[11] + a.data[27],
        a.data[12] + a.data[28],
        a.data[13] + a.data[29],
        a.data[14] + a.data[30],
        a.data[15] + a.data[31]});
    splitcomplex4 q = splitcomplex4_make((const double[]){
        a.data[0] - a.data[16],
        a.data[1] - a.data[17],
        a.data[2] - a.data[18],
        a.data[3] - a.data[19],
        a.data[4] - a.data[20],
        a.data[5] - a.data[21],
        a.data[6] - a.data[22],
        a.data[7] - a.data[23],
        a.data[8] - a.data[24],
        a.data[9] - a.data[25],
        a.data[10] - a.data[26],
        a.data[11] - a.data[27],
        a.data[12] - a.data[28],
        a.data[13] - a.data[29],
        a.data[14] - a.data[30],
        a.data[15] - a.data[31]});
    p = f(p);
    q = f(q);
    return splitcomplex5_make((const double[]){
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
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (p.data[2] - q.data[2]) / 2,
        (p.data[3] - q.data[3]) / 2,
        (p.data[4] - q.data[4]) / 2,
        (p.data[5] - q.data[5]) / 2,
        (p.data[6] - q.data[6]) / 2,
        (p.data[7] - q.data[7]) / 2,
        (p.data[8] - q.data[8]) / 2,
        (p.data[9] - q.data[9]) / 2,
        (p.data[10] - q.data[10]) / 2,
        (p.data[11] - q.data[11]) / 2,
        (p.data[12] - q.data[12]) / 2,
        (p.data[13] - q.data[13]) / 2,
        (p.data[14] - q.data[14]) / 2,
        (p.data[15] - q.data[15]) / 2});
}

static inline splitcomplex5
splitcomplex5_from_re(double a)
{
    splitcomplex5 result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
splitcomplex5_to_re(splitcomplex5 a)
{
    return a.data[0];
}

static inline splitcomplex5
splitcomplex5_zero(void)
{
    return vec32_zero();
}

static inline splitcomplex5
splitcomplex5_one(void)
{
    return splitcomplex5_from_re(1);
}

static inline int
splitcomplex5_is_real(splitcomplex5 a)
{
    return vec32_is_real(a);
}

static inline splitcomplex5
splitcomplex5_pos(splitcomplex5 a)
{
    return vec32_pos(a);
}

static inline splitcomplex5
splitcomplex5_neg(splitcomplex5 a)
{
    return vec32_neg(a);
}

static inline splitcomplex5
splitcomplex5_re(splitcomplex5 a)
{
    return splitcomplex5_from_re(a.data[0]);
}

static inline splitcomplex5
splitcomplex5_im(splitcomplex5 a)
{
    return vec32_sub(a, splitcomplex5_re(a));
}

static inline splitcomplex5
splitcomplex5_conjugate(splitcomplex5 a)
{
    return splitcomplex5_make((const double[]){a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15], -a.data[16], -a.data[17], -a.data[18], -a.data[19], -a.data[20], -a.data[21], -a.data[22], -a.data[23], -a.data[24], -a.data[25], -a.data[26], -a.data[27], -a.data[28], -a.data[29], -a.data[30], -a.data[31]});
}

static inline splitcomplex5
splitcomplex5_add(splitcomplex5 a, splitcomplex5 b)
{
    return vec32_add(a, b);
}

static inline splitcomplex5
splitcomplex5_sub(splitcomplex5 a, splitcomplex5 b)
{
    return vec32_sub(a, b);
}

static inline int
splitcomplex5_equal(splitcomplex5 a, splitcomplex5 b)
{
    return vec32_all_eq(a, b);
}

static inline int
splitcomplex5_not_equal(splitcomplex5 a, splitcomplex5 b)
{
    return !splitcomplex5_equal(a, b);
}

static inline splitcomplex5
splitcomplex5_mul(splitcomplex5 a, splitcomplex5 b)
{
    if (splitcomplex5_is_real(b)) {
        return vec32_scale(a, b.data[0]);
    }
    if (splitcomplex5_is_real(a)) {
        return vec32_scale(b, a.data[0]);
    }
    splitcomplex4 u = splitcomplex4_make(a.data);
    splitcomplex4 v = splitcomplex4_make(a.data + 16);
    splitcomplex4 s = splitcomplex4_make(b.data);
    splitcomplex4 t = splitcomplex4_make(b.data + 16);
    splitcomplex4 lo = splitcomplex4_add(splitcomplex4_mul(u, s), splitcomplex4_mul(v, t));
    splitcomplex4 hi = splitcomplex4_add(splitcomplex4_mul(u, t), splitcomplex4_mul(v, s));
    splitcomplex5 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_SPLITCOMPLEX5_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

/* Recursive norm: (u + v*e)*(u - v*e) = u*u - v*v. */
static inline splitcomplex4
splitcomplex5_norm2(splitcomplex5 a)
{
    splitcomplex4 u = splitcomplex4_make(a.data);
    splitcomplex4 v = splitcomplex4_make(a.data + 16);
    return splitcomplex4_sub(splitcomplex4_mul(u, u), splitcomplex4_mul(v, v));
}

static inline splitcomplex3
splitcomplex5_norm4(splitcomplex5 a)
{
    return splitcomplex4_norm2(splitcomplex5_norm2(a));
}

static inline splitcomplex2
splitcomplex5_norm8(splitcomplex5 a)
{
    return splitcomplex3_norm2(splitcomplex5_norm4(a));
}

static inline splitcomplex1
splitcomplex5_norm16(splitcomplex5 a)
{
    return splitcomplex2_norm2(splitcomplex5_norm8(a));
}

static inline double
splitcomplex5_norm32(splitcomplex5 a)
{
    return splitcomplex1_norm2(splitcomplex5_norm16(a));
}

static inline complex
splitcomplex5_norm(splitcomplex5 a)
{
    return complex_root_real(splitcomplex5_norm32(a), 32);
}

static inline splitcomplex5
splitcomplex5_inv(splitcomplex5 a)
{
    splitcomplex4 ninv = splitcomplex4_inv(splitcomplex5_norm2(a));
    splitcomplex4 u = splitcomplex4_make(a.data);
    splitcomplex4 v = splitcomplex4_make(a.data + 16);
    splitcomplex4 lo = splitcomplex4_mul(u, ninv);
    splitcomplex4 hi = splitcomplex4_mul(splitcomplex4_neg(v), ninv);
    splitcomplex5 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_SPLITCOMPLEX5_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

static inline splitcomplex5
splitcomplex5_div(splitcomplex5 a, splitcomplex5 b)
{
    if (splitcomplex5_is_real(b)) {
        return vec32_div_scalar(a, b.data[0]);
    }
    return splitcomplex5_mul(a, splitcomplex5_inv(b));
}

static inline splitcomplex5
splitcomplex5_square(splitcomplex5 a)
{
    return splitcomplex5_mul(a, a);
}

static inline splitcomplex5
splitcomplex5_exp(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_exp);
}

static inline splitcomplex5
splitcomplex5_log(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_log);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline splitcomplex5
splitcomplex5_ipow(splitcomplex5 a, int64_t b)
{
    splitcomplex5 result = splitcomplex5_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return splitcomplex5_mul(a, a);
    }
    if (b == -1) {
        return splitcomplex5_inv(a);
    }
    if (b < 0) {
        a = splitcomplex5_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = splitcomplex5_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = splitcomplex5_mul(a, a);
        }
    }
    return result;
}

static inline splitcomplex5
splitcomplex5_rpow(splitcomplex5 a, splitcomplex5 b)
{
    splitcomplex4 ap = splitcomplex4_make((const double[]){
        a.data[0] + a.data[16],
        a.data[1] + a.data[17],
        a.data[2] + a.data[18],
        a.data[3] + a.data[19],
        a.data[4] + a.data[20],
        a.data[5] + a.data[21],
        a.data[6] + a.data[22],
        a.data[7] + a.data[23],
        a.data[8] + a.data[24],
        a.data[9] + a.data[25],
        a.data[10] + a.data[26],
        a.data[11] + a.data[27],
        a.data[12] + a.data[28],
        a.data[13] + a.data[29],
        a.data[14] + a.data[30],
        a.data[15] + a.data[31]});
    splitcomplex4 aq = splitcomplex4_make((const double[]){
        a.data[0] - a.data[16],
        a.data[1] - a.data[17],
        a.data[2] - a.data[18],
        a.data[3] - a.data[19],
        a.data[4] - a.data[20],
        a.data[5] - a.data[21],
        a.data[6] - a.data[22],
        a.data[7] - a.data[23],
        a.data[8] - a.data[24],
        a.data[9] - a.data[25],
        a.data[10] - a.data[26],
        a.data[11] - a.data[27],
        a.data[12] - a.data[28],
        a.data[13] - a.data[29],
        a.data[14] - a.data[30],
        a.data[15] - a.data[31]});
    splitcomplex4 bp = splitcomplex4_make((const double[]){
        b.data[0] + b.data[16],
        b.data[1] + b.data[17],
        b.data[2] + b.data[18],
        b.data[3] + b.data[19],
        b.data[4] + b.data[20],
        b.data[5] + b.data[21],
        b.data[6] + b.data[22],
        b.data[7] + b.data[23],
        b.data[8] + b.data[24],
        b.data[9] + b.data[25],
        b.data[10] + b.data[26],
        b.data[11] + b.data[27],
        b.data[12] + b.data[28],
        b.data[13] + b.data[29],
        b.data[14] + b.data[30],
        b.data[15] + b.data[31]});
    splitcomplex4 bq = splitcomplex4_make((const double[]){
        b.data[0] - b.data[16],
        b.data[1] - b.data[17],
        b.data[2] - b.data[18],
        b.data[3] - b.data[19],
        b.data[4] - b.data[20],
        b.data[5] - b.data[21],
        b.data[6] - b.data[22],
        b.data[7] - b.data[23],
        b.data[8] - b.data[24],
        b.data[9] - b.data[25],
        b.data[10] - b.data[26],
        b.data[11] - b.data[27],
        b.data[12] - b.data[28],
        b.data[13] - b.data[29],
        b.data[14] - b.data[30],
        b.data[15] - b.data[31]});
    splitcomplex4 p = splitcomplex4_rpow(ap, bp);
    splitcomplex4 q = splitcomplex4_rpow(aq, bq);
    return splitcomplex5_make((const double[]){
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
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (p.data[2] - q.data[2]) / 2,
        (p.data[3] - q.data[3]) / 2,
        (p.data[4] - q.data[4]) / 2,
        (p.data[5] - q.data[5]) / 2,
        (p.data[6] - q.data[6]) / 2,
        (p.data[7] - q.data[7]) / 2,
        (p.data[8] - q.data[8]) / 2,
        (p.data[9] - q.data[9]) / 2,
        (p.data[10] - q.data[10]) / 2,
        (p.data[11] - q.data[11]) / 2,
        (p.data[12] - q.data[12]) / 2,
        (p.data[13] - q.data[13]) / 2,
        (p.data[14] - q.data[14]) / 2,
        (p.data[15] - q.data[15]) / 2});
}

static inline splitcomplex5
splitcomplex5_pow(splitcomplex5 a, splitcomplex5 b)
{
    return splitcomplex5_rpow(a, b);
}

static inline splitcomplex5
splitcomplex5_lpow(splitcomplex5 a, splitcomplex5 b)
{
    return splitcomplex5_rpow(a, b);
}

static inline splitcomplex5
splitcomplex5_sqrt(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_sqrt);
}

static inline splitcomplex5
splitcomplex5_cos(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_cos);
}

static inline splitcomplex5
splitcomplex5_sin(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_sin);
}

static inline splitcomplex5
splitcomplex5_tan(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_tan);
}

static inline splitcomplex5
splitcomplex5_acos(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_acos);
}

static inline splitcomplex5
splitcomplex5_asin(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_asin);
}

static inline splitcomplex5
splitcomplex5_atan(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_atan);
}

static inline splitcomplex5
splitcomplex5_cosh(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_cosh);
}

static inline splitcomplex5
splitcomplex5_sinh(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_sinh);
}

static inline splitcomplex5
splitcomplex5_tanh(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_tanh);
}

static inline splitcomplex5
splitcomplex5_acosh(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_acosh);
}

static inline splitcomplex5
splitcomplex5_asinh(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_asinh);
}

static inline splitcomplex5
splitcomplex5_atanh(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_atanh);
}

static inline splitcomplex5
splitcomplex5_expm1(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_expm1);
}

static inline splitcomplex5
splitcomplex5_log1p(splitcomplex5 a)
{
    return splitcomplex5_apply(a, splitcomplex4_log1p);
}

static inline splitcomplex5
splitcomplex5_radd(double a, splitcomplex5 b)
{
    return splitcomplex5_add(splitcomplex5_from_re(a), b);
}

static inline splitcomplex5
splitcomplex5_ladd(splitcomplex5 a, double b)
{
    return splitcomplex5_add(a, splitcomplex5_from_re(b));
}

static inline splitcomplex5
splitcomplex5_rsub(double a, splitcomplex5 b)
{
    return splitcomplex5_sub(splitcomplex5_from_re(a), b);
}

static inline splitcomplex5
splitcomplex5_lsub(splitcomplex5 a, double b)
{
    return splitcomplex5_sub(a, splitcomplex5_from_re(b));
}

static inline splitcomplex5
splitcomplex5_rmul(double a, splitcomplex5 b)
{
    return vec32_scale(b, a);
}

static inline splitcomplex5
splitcomplex5_lmul(splitcomplex5 a, double b)
{
    return vec32_scale(a, b);
}

static inline splitcomplex5
splitcomplex5_rdiv(double a, splitcomplex5 b)
{
    return splitcomplex5_div(splitcomplex5_from_re(a), b);
}

static inline splitcomplex5
splitcomplex5_ldiv(splitcomplex5 a, double b)
{
    return vec32_div_scalar(a, b);
}

static inline splitcomplex5
splitcomplex5_rpows(double a, splitcomplex5 b)
{
    return splitcomplex5_rpow(splitcomplex5_from_re(a), b);
}

static inline splitcomplex5
splitcomplex5_lpows(splitcomplex5 a, double b)
{
    return splitcomplex5_rpow(a, splitcomplex5_from_re(b));
}

static inline splitcomplex5
splitcomplex5_mix(splitcomplex5 a, splitcomplex5 b, double t)
{
    splitcomplex5 delta = splitcomplex5_mul(splitcomplex5_inv(a), b);
    return splitcomplex5_mul(a, splitcomplex5_lpows(delta, t));
}

static inline void
splitcomplex5_basis(splitcomplex5 result[32])
{
    int i, j;
    for (i = 0; i < ALGPY_SPLITCOMPLEX5_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX5_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
splitcomplex5_gens(splitcomplex5 result[5])
{
    result[0] = splitcomplex5_make((const double[ALGPY_SPLITCOMPLEX5_N]){[1] = 1});
    result[1] = splitcomplex5_make((const double[ALGPY_SPLITCOMPLEX5_N]){[2] = 1});
    result[2] = splitcomplex5_make((const double[ALGPY_SPLITCOMPLEX5_N]){[4] = 1});
    result[3] = splitcomplex5_make((const double[ALGPY_SPLITCOMPLEX5_N]){[8] = 1});
    result[4] = splitcomplex5_make((const double[ALGPY_SPLITCOMPLEX5_N]){[16] = 1});
}

static inline void
splitcomplex5_coeffs(splitcomplex5 result[ALGPY_SPLITCOMPLEX5_N*ALGPY_SPLITCOMPLEX5_N])
{
    splitcomplex5 basis[ALGPY_SPLITCOMPLEX5_N];
    int i, j;
    splitcomplex5_basis(basis);
    for (i = 0; i < ALGPY_SPLITCOMPLEX5_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX5_N; j++) {
            result[i*ALGPY_SPLITCOMPLEX5_N+j] = splitcomplex5_mul(basis[i], basis[j]);
        }
    }
}

static inline void
splitcomplex5_minpoly(splitcomplex5 a, double result[33])
{
    splitcomplex4 p = splitcomplex4_make((const double[]){
        a.data[0] + a.data[16],
        a.data[1] + a.data[17],
        a.data[2] + a.data[18],
        a.data[3] + a.data[19],
        a.data[4] + a.data[20],
        a.data[5] + a.data[21],
        a.data[6] + a.data[22],
        a.data[7] + a.data[23],
        a.data[8] + a.data[24],
        a.data[9] + a.data[25],
        a.data[10] + a.data[26],
        a.data[11] + a.data[27],
        a.data[12] + a.data[28],
        a.data[13] + a.data[29],
        a.data[14] + a.data[30],
        a.data[15] + a.data[31]});
    splitcomplex4 q = splitcomplex4_make((const double[]){
        a.data[0] - a.data[16],
        a.data[1] - a.data[17],
        a.data[2] - a.data[18],
        a.data[3] - a.data[19],
        a.data[4] - a.data[20],
        a.data[5] - a.data[21],
        a.data[6] - a.data[22],
        a.data[7] - a.data[23],
        a.data[8] - a.data[24],
        a.data[9] - a.data[25],
        a.data[10] - a.data[26],
        a.data[11] - a.data[27],
        a.data[12] - a.data[28],
        a.data[13] - a.data[29],
        a.data[14] - a.data[30],
        a.data[15] - a.data[31]});
    double left[17], right[17];
    int i, j;
    splitcomplex4_minpoly(p, left);
    splitcomplex4_minpoly(q, right);
    for (i = 0; i <= 32; i++) result[i] = 0;
    for (i = 0; i <= 16; i++) {
        for (j = 0; j <= 16; j++) result[i+j] += left[i]*right[j];
    }
}

static inline void
splitcomplex5_eivals(splitcomplex5 a, complex result[32])
{
    splitcomplex4 p = splitcomplex4_make((const double[]){
        a.data[0] + a.data[16],
        a.data[1] + a.data[17],
        a.data[2] + a.data[18],
        a.data[3] + a.data[19],
        a.data[4] + a.data[20],
        a.data[5] + a.data[21],
        a.data[6] + a.data[22],
        a.data[7] + a.data[23],
        a.data[8] + a.data[24],
        a.data[9] + a.data[25],
        a.data[10] + a.data[26],
        a.data[11] + a.data[27],
        a.data[12] + a.data[28],
        a.data[13] + a.data[29],
        a.data[14] + a.data[30],
        a.data[15] + a.data[31]});
    splitcomplex4 q = splitcomplex4_make((const double[]){
        a.data[0] - a.data[16],
        a.data[1] - a.data[17],
        a.data[2] - a.data[18],
        a.data[3] - a.data[19],
        a.data[4] - a.data[20],
        a.data[5] - a.data[21],
        a.data[6] - a.data[22],
        a.data[7] - a.data[23],
        a.data[8] - a.data[24],
        a.data[9] - a.data[25],
        a.data[10] - a.data[26],
        a.data[11] - a.data[27],
        a.data[12] - a.data[28],
        a.data[13] - a.data[29],
        a.data[14] - a.data[30],
        a.data[15] - a.data[31]});
    complex left[16], right[16];
    splitcomplex4_eivals(p, left);
    splitcomplex4_eivals(q, right);
    /* Preserve the sign-bit ordering of the generators. */
    result[0] = right[0];
    result[1] = right[1];
    result[2] = right[2];
    result[3] = right[3];
    result[4] = right[4];
    result[5] = right[5];
    result[6] = right[6];
    result[7] = right[7];
    result[8] = right[8];
    result[9] = right[9];
    result[10] = right[10];
    result[11] = right[11];
    result[12] = right[12];
    result[13] = right[13];
    result[14] = right[14];
    result[15] = right[15];
    result[16] = left[0];
    result[17] = left[1];
    result[18] = left[2];
    result[19] = left[3];
    result[20] = left[4];
    result[21] = left[5];
    result[22] = left[6];
    result[23] = left[7];
    result[24] = left[8];
    result[25] = left[9];
    result[26] = left[10];
    result[27] = left[11];
    result[28] = left[12];
    result[29] = left[13];
    result[30] = left[14];
    result[31] = left[15];
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
splitcomplex5_lmat(splitcomplex5 a, double result[ALGPY_SPLITCOMPLEX5_N*ALGPY_SPLITCOMPLEX5_N])
{
    splitcomplex5 basis = {{0}};
    splitcomplex5 product;
    int row, col;
    for (col = 0; col < ALGPY_SPLITCOMPLEX5_N; col++) {
        basis.data[col] = 1;
        product = splitcomplex5_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_SPLITCOMPLEX5_N; row++) {
            result[row * ALGPY_SPLITCOMPLEX5_N + col] = product.data[row];
        }
    }
}

static inline void
splitcomplex5_rmat(splitcomplex5 a, double result[ALGPY_SPLITCOMPLEX5_N*ALGPY_SPLITCOMPLEX5_N])
{
    splitcomplex5 basis = {{0}};
    splitcomplex5 product;
    int row, col;
    for (row = 0; row < ALGPY_SPLITCOMPLEX5_N; row++) {
        basis.data[row] = 1;
        product = splitcomplex5_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_SPLITCOMPLEX5_N; col++) {
            result[row * ALGPY_SPLITCOMPLEX5_N + col] = product.data[col];
        }
    }
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline splitcomplex5
splitcomplex5_grade0(splitcomplex5 a)
{
    const double *x = a.data;
    return splitcomplex5_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline splitcomplex5
splitcomplex5_grade1(splitcomplex5 a)
{
    const double *x = a.data;
    return splitcomplex5_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline splitcomplex5
splitcomplex5_grade2(splitcomplex5 a)
{
    const double *x = a.data;
    return splitcomplex5_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline splitcomplex5
splitcomplex5_grade3(splitcomplex5 a)
{
    const double *x = a.data;
    return splitcomplex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline splitcomplex5
splitcomplex5_grade4(splitcomplex5 a)
{
    const double *x = a.data;
    return splitcomplex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline splitcomplex5
splitcomplex5_grade5(splitcomplex5 a)
{
    const double *x = a.data;
    return splitcomplex5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

#endif /* ALGPY_SPLITCOMPLEX5_H */
