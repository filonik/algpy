#ifndef ALGPY_SPLITCOMPLEX4_H
#define ALGPY_SPLITCOMPLEX4_H

#include <stdint.h>
#include "../../vectorspaces/vec16.h"
#include "splitcomplex3.h"

#define ALGPY_SPLITCOMPLEX4_N ALGPY_VEC16_N

typedef vec16 splitcomplex4;

static inline splitcomplex4
splitcomplex4_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10, double x11, double x12, double x13, double x14, double x15)
{
    splitcomplex4 result = {{x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15}};
    return result;
}

/* For a = u + v*e with e*e = +1, apply f to u+v and u-v. */
static inline splitcomplex4
splitcomplex4_apply(splitcomplex4 a, splitcomplex3 (*f)(splitcomplex3))
{
    splitcomplex3 p = splitcomplex3_make(
        a.data[0] + a.data[8],
        a.data[1] + a.data[9],
        a.data[2] + a.data[10],
        a.data[3] + a.data[11],
        a.data[4] + a.data[12],
        a.data[5] + a.data[13],
        a.data[6] + a.data[14],
        a.data[7] + a.data[15]);
    splitcomplex3 q = splitcomplex3_make(
        a.data[0] - a.data[8],
        a.data[1] - a.data[9],
        a.data[2] - a.data[10],
        a.data[3] - a.data[11],
        a.data[4] - a.data[12],
        a.data[5] - a.data[13],
        a.data[6] - a.data[14],
        a.data[7] - a.data[15]);
    p = f(p);
    q = f(q);
    return splitcomplex4_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[2] + q.data[2]) / 2,
        (p.data[3] + q.data[3]) / 2,
        (p.data[4] + q.data[4]) / 2,
        (p.data[5] + q.data[5]) / 2,
        (p.data[6] + q.data[6]) / 2,
        (p.data[7] + q.data[7]) / 2,
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (p.data[2] - q.data[2]) / 2,
        (p.data[3] - q.data[3]) / 2,
        (p.data[4] - q.data[4]) / 2,
        (p.data[5] - q.data[5]) / 2,
        (p.data[6] - q.data[6]) / 2,
        (p.data[7] - q.data[7]) / 2);
}

static inline splitcomplex4
splitcomplex4_from_re(double a)
{
    return splitcomplex4_make(a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
splitcomplex4_to_re(splitcomplex4 a)
{
    return a.data[0];
}

static inline splitcomplex4
splitcomplex4_zero(void)
{
    return vec16_zero();
}

static inline splitcomplex4
splitcomplex4_one(void)
{
    return splitcomplex4_from_re(1);
}

static inline int
splitcomplex4_is_real(splitcomplex4 a)
{
    return vec16_is_real(a);
}

static inline splitcomplex4
splitcomplex4_pos(splitcomplex4 a)
{
    return vec16_pos(a);
}

static inline splitcomplex4
splitcomplex4_neg(splitcomplex4 a)
{
    return vec16_neg(a);
}

static inline splitcomplex4
splitcomplex4_re(splitcomplex4 a)
{
    return splitcomplex4_from_re(a.data[0]);
}

static inline splitcomplex4
splitcomplex4_im(splitcomplex4 a)
{
    return vec16_sub(a, splitcomplex4_re(a));
}

static inline splitcomplex4
splitcomplex4_conjugate(splitcomplex4 a)
{
    return splitcomplex4_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], -a.data[8], -a.data[9], -a.data[10], -a.data[11], -a.data[12], -a.data[13], -a.data[14], -a.data[15]);
}

static inline splitcomplex4
splitcomplex4_add(splitcomplex4 a, splitcomplex4 b)
{
    return vec16_add(a, b);
}

static inline splitcomplex4
splitcomplex4_sub(splitcomplex4 a, splitcomplex4 b)
{
    return vec16_sub(a, b);
}

static inline int
splitcomplex4_equal(splitcomplex4 a, splitcomplex4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
splitcomplex4_not_equal(splitcomplex4 a, splitcomplex4 b)
{
    return !splitcomplex4_equal(a, b);
}

static inline splitcomplex4
splitcomplex4_mul(splitcomplex4 a, splitcomplex4 b)
{
    if (splitcomplex4_is_real(b)) {
        return vec16_scale(a, b.data[0]);
    }
    if (splitcomplex4_is_real(a)) {
        return vec16_scale(b, a.data[0]);
    }
    splitcomplex3 u = splitcomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    splitcomplex3 v = splitcomplex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    splitcomplex3 s = splitcomplex3_make(b.data[0], b.data[1], b.data[2], b.data[3], b.data[4], b.data[5], b.data[6], b.data[7]);
    splitcomplex3 t = splitcomplex3_make(b.data[8], b.data[9], b.data[10], b.data[11], b.data[12], b.data[13], b.data[14], b.data[15]);
    splitcomplex3 lo = splitcomplex3_add(splitcomplex3_mul(u, s), splitcomplex3_mul(v, t));
    splitcomplex3 hi = splitcomplex3_add(splitcomplex3_mul(u, t), splitcomplex3_mul(v, s));
    return splitcomplex4_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], lo.data[4], lo.data[5], lo.data[6], lo.data[7], hi.data[0], hi.data[1], hi.data[2], hi.data[3], hi.data[4], hi.data[5], hi.data[6], hi.data[7]);
}

/* Recursive norm: (u + v*e)*(u - v*e) = u*u - v*v. */
static inline splitcomplex3
splitcomplex4_norm2(splitcomplex4 a)
{
    splitcomplex3 u = splitcomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    splitcomplex3 v = splitcomplex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    return splitcomplex3_sub(splitcomplex3_mul(u, u), splitcomplex3_mul(v, v));
}

static inline splitcomplex2
splitcomplex4_norm4(splitcomplex4 a)
{
    return splitcomplex3_norm2(splitcomplex4_norm2(a));
}

static inline splitcomplex1
splitcomplex4_norm8(splitcomplex4 a)
{
    return splitcomplex2_norm2(splitcomplex4_norm4(a));
}

static inline double
splitcomplex4_norm16(splitcomplex4 a)
{
    return splitcomplex1_norm2(splitcomplex4_norm8(a));
}

static inline complex
splitcomplex4_norm(splitcomplex4 a)
{
    return complex_root_real(splitcomplex4_norm16(a), 16);
}

static inline splitcomplex4
splitcomplex4_inv(splitcomplex4 a)
{
    splitcomplex3 ninv = splitcomplex3_inv(splitcomplex4_norm2(a));
    splitcomplex3 u = splitcomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    splitcomplex3 v = splitcomplex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    splitcomplex3 lo = splitcomplex3_mul(u, ninv);
    splitcomplex3 hi = splitcomplex3_mul(splitcomplex3_neg(v), ninv);
    return splitcomplex4_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], lo.data[4], lo.data[5], lo.data[6], lo.data[7], hi.data[0], hi.data[1], hi.data[2], hi.data[3], hi.data[4], hi.data[5], hi.data[6], hi.data[7]);
}

static inline splitcomplex4
splitcomplex4_div(splitcomplex4 a, splitcomplex4 b)
{
    if (splitcomplex4_is_real(b)) {
        return vec16_div_scalar(a, b.data[0]);
    }
    return splitcomplex4_mul(a, splitcomplex4_inv(b));
}

static inline splitcomplex4
splitcomplex4_square(splitcomplex4 a)
{
    return splitcomplex4_mul(a, a);
}

static inline splitcomplex4
splitcomplex4_exp(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_exp);
}

static inline splitcomplex4
splitcomplex4_log(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_log);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline splitcomplex4
splitcomplex4_ipow(splitcomplex4 a, int64_t b)
{
    splitcomplex4 result = splitcomplex4_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return splitcomplex4_mul(a, a);
    }
    if (b == -1) {
        return splitcomplex4_inv(a);
    }
    if (b < 0) {
        a = splitcomplex4_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = splitcomplex4_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = splitcomplex4_mul(a, a);
        }
    }
    return result;
}

static inline splitcomplex4
splitcomplex4_rpow(splitcomplex4 a, splitcomplex4 b)
{
    splitcomplex3 ap = splitcomplex3_make(
        a.data[0] + a.data[8],
        a.data[1] + a.data[9],
        a.data[2] + a.data[10],
        a.data[3] + a.data[11],
        a.data[4] + a.data[12],
        a.data[5] + a.data[13],
        a.data[6] + a.data[14],
        a.data[7] + a.data[15]);
    splitcomplex3 aq = splitcomplex3_make(
        a.data[0] - a.data[8],
        a.data[1] - a.data[9],
        a.data[2] - a.data[10],
        a.data[3] - a.data[11],
        a.data[4] - a.data[12],
        a.data[5] - a.data[13],
        a.data[6] - a.data[14],
        a.data[7] - a.data[15]);
    splitcomplex3 bp = splitcomplex3_make(
        b.data[0] + b.data[8],
        b.data[1] + b.data[9],
        b.data[2] + b.data[10],
        b.data[3] + b.data[11],
        b.data[4] + b.data[12],
        b.data[5] + b.data[13],
        b.data[6] + b.data[14],
        b.data[7] + b.data[15]);
    splitcomplex3 bq = splitcomplex3_make(
        b.data[0] - b.data[8],
        b.data[1] - b.data[9],
        b.data[2] - b.data[10],
        b.data[3] - b.data[11],
        b.data[4] - b.data[12],
        b.data[5] - b.data[13],
        b.data[6] - b.data[14],
        b.data[7] - b.data[15]);
    splitcomplex3 p = splitcomplex3_rpow(ap, bp);
    splitcomplex3 q = splitcomplex3_rpow(aq, bq);
    return splitcomplex4_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[2] + q.data[2]) / 2,
        (p.data[3] + q.data[3]) / 2,
        (p.data[4] + q.data[4]) / 2,
        (p.data[5] + q.data[5]) / 2,
        (p.data[6] + q.data[6]) / 2,
        (p.data[7] + q.data[7]) / 2,
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (p.data[2] - q.data[2]) / 2,
        (p.data[3] - q.data[3]) / 2,
        (p.data[4] - q.data[4]) / 2,
        (p.data[5] - q.data[5]) / 2,
        (p.data[6] - q.data[6]) / 2,
        (p.data[7] - q.data[7]) / 2);
}

static inline splitcomplex4
splitcomplex4_pow(splitcomplex4 a, splitcomplex4 b)
{
    return splitcomplex4_rpow(a, b);
}

static inline splitcomplex4
splitcomplex4_lpow(splitcomplex4 a, splitcomplex4 b)
{
    return splitcomplex4_rpow(a, b);
}

static inline splitcomplex4
splitcomplex4_sqrt(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_sqrt);
}

static inline splitcomplex4
splitcomplex4_cos(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_cos);
}

static inline splitcomplex4
splitcomplex4_sin(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_sin);
}

static inline splitcomplex4
splitcomplex4_tan(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_tan);
}

static inline splitcomplex4
splitcomplex4_acos(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_acos);
}

static inline splitcomplex4
splitcomplex4_asin(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_asin);
}

static inline splitcomplex4
splitcomplex4_atan(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_atan);
}

static inline splitcomplex4
splitcomplex4_cosh(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_cosh);
}

static inline splitcomplex4
splitcomplex4_sinh(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_sinh);
}

static inline splitcomplex4
splitcomplex4_tanh(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_tanh);
}

static inline splitcomplex4
splitcomplex4_acosh(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_acosh);
}

static inline splitcomplex4
splitcomplex4_asinh(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_asinh);
}

static inline splitcomplex4
splitcomplex4_atanh(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_atanh);
}

static inline splitcomplex4
splitcomplex4_expm1(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_expm1);
}

static inline splitcomplex4
splitcomplex4_log1p(splitcomplex4 a)
{
    return splitcomplex4_apply(a, splitcomplex3_log1p);
}

static inline splitcomplex4
splitcomplex4_radd(double a, splitcomplex4 b)
{
    return splitcomplex4_add(splitcomplex4_from_re(a), b);
}

static inline splitcomplex4
splitcomplex4_ladd(splitcomplex4 a, double b)
{
    return splitcomplex4_add(a, splitcomplex4_from_re(b));
}

static inline splitcomplex4
splitcomplex4_rsub(double a, splitcomplex4 b)
{
    return splitcomplex4_sub(splitcomplex4_from_re(a), b);
}

static inline splitcomplex4
splitcomplex4_lsub(splitcomplex4 a, double b)
{
    return splitcomplex4_sub(a, splitcomplex4_from_re(b));
}

static inline splitcomplex4
splitcomplex4_rmul(double a, splitcomplex4 b)
{
    return vec16_scale(b, a);
}

static inline splitcomplex4
splitcomplex4_lmul(splitcomplex4 a, double b)
{
    return vec16_scale(a, b);
}

static inline splitcomplex4
splitcomplex4_rdiv(double a, splitcomplex4 b)
{
    return splitcomplex4_div(splitcomplex4_from_re(a), b);
}

static inline splitcomplex4
splitcomplex4_ldiv(splitcomplex4 a, double b)
{
    return vec16_div_scalar(a, b);
}

static inline splitcomplex4
splitcomplex4_rpows(double a, splitcomplex4 b)
{
    return splitcomplex4_rpow(splitcomplex4_from_re(a), b);
}

static inline splitcomplex4
splitcomplex4_lpows(splitcomplex4 a, double b)
{
    return splitcomplex4_rpow(a, splitcomplex4_from_re(b));
}

static inline splitcomplex4
splitcomplex4_mix(splitcomplex4 a, splitcomplex4 b, double t)
{
    splitcomplex4 delta = splitcomplex4_mul(splitcomplex4_inv(a), b);
    return splitcomplex4_mul(a, splitcomplex4_lpows(delta, t));
}

static inline void
splitcomplex4_basis(splitcomplex4 result[16])
{
    int i, j;
    for (i = 0; i < ALGPY_SPLITCOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
splitcomplex4_coeffs(splitcomplex4 result[ALGPY_SPLITCOMPLEX4_N*ALGPY_SPLITCOMPLEX4_N])
{
    splitcomplex4 basis[ALGPY_SPLITCOMPLEX4_N];
    int i, j;
    splitcomplex4_basis(basis);
    for (i = 0; i < ALGPY_SPLITCOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX4_N; j++) {
            result[i*ALGPY_SPLITCOMPLEX4_N+j] = splitcomplex4_mul(basis[i], basis[j]);
        }
    }
}

static inline void
splitcomplex4_minpoly(splitcomplex4 a, double result[17])
{
    splitcomplex3 p = splitcomplex3_make(
        a.data[0] + a.data[8],
        a.data[1] + a.data[9],
        a.data[2] + a.data[10],
        a.data[3] + a.data[11],
        a.data[4] + a.data[12],
        a.data[5] + a.data[13],
        a.data[6] + a.data[14],
        a.data[7] + a.data[15]);
    splitcomplex3 q = splitcomplex3_make(
        a.data[0] - a.data[8],
        a.data[1] - a.data[9],
        a.data[2] - a.data[10],
        a.data[3] - a.data[11],
        a.data[4] - a.data[12],
        a.data[5] - a.data[13],
        a.data[6] - a.data[14],
        a.data[7] - a.data[15]);
    double left[9], right[9];
    int i, j;
    splitcomplex3_minpoly(p, left);
    splitcomplex3_minpoly(q, right);
    for (i = 0; i <= 16; i++) {
        result[i] = 0;
    }
    for (i = 0; i <= 8; i++) {
        for (j = 0; j <= 8; j++) {
            result[i+j] += left[i]*right[j];
        }
    }
}

static inline void
splitcomplex4_eivals(splitcomplex4 a, complex result[16])
{
    splitcomplex3 p = splitcomplex3_make(
        a.data[0] + a.data[8],
        a.data[1] + a.data[9],
        a.data[2] + a.data[10],
        a.data[3] + a.data[11],
        a.data[4] + a.data[12],
        a.data[5] + a.data[13],
        a.data[6] + a.data[14],
        a.data[7] + a.data[15]);
    splitcomplex3 q = splitcomplex3_make(
        a.data[0] - a.data[8],
        a.data[1] - a.data[9],
        a.data[2] - a.data[10],
        a.data[3] - a.data[11],
        a.data[4] - a.data[12],
        a.data[5] - a.data[13],
        a.data[6] - a.data[14],
        a.data[7] - a.data[15]);
    complex left[8], right[8];
    splitcomplex3_eivals(p, left);
    splitcomplex3_eivals(q, right);
    /* Preserve the sign-bit ordering of the generators. */
    result[0] = right[0];
    result[1] = right[1];
    result[2] = right[2];
    result[3] = right[3];
    result[4] = right[4];
    result[5] = right[5];
    result[6] = right[6];
    result[7] = right[7];
    result[8] = left[0];
    result[9] = left[1];
    result[10] = left[2];
    result[11] = left[3];
    result[12] = left[4];
    result[13] = left[5];
    result[14] = left[6];
    result[15] = left[7];
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
splitcomplex4_lmat(splitcomplex4 a, double result[ALGPY_SPLITCOMPLEX4_N*ALGPY_SPLITCOMPLEX4_N])
{
    splitcomplex4 basis = {{0}};
    splitcomplex4 product;
    int row, col;
    for (col = 0; col < ALGPY_SPLITCOMPLEX4_N; col++) {
        basis.data[col] = 1;
        product = splitcomplex4_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_SPLITCOMPLEX4_N; row++) {
            result[row * ALGPY_SPLITCOMPLEX4_N + col] = product.data[row];
        }
    }
}

static inline void
splitcomplex4_rmat(splitcomplex4 a, double result[ALGPY_SPLITCOMPLEX4_N*ALGPY_SPLITCOMPLEX4_N])
{
    splitcomplex4 basis = {{0}};
    splitcomplex4 product;
    int row, col;
    for (row = 0; row < ALGPY_SPLITCOMPLEX4_N; row++) {
        basis.data[row] = 1;
        product = splitcomplex4_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_SPLITCOMPLEX4_N; col++) {
            result[row * ALGPY_SPLITCOMPLEX4_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_SPLITCOMPLEX4_H */
