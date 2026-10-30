#ifndef ALGPY_COMPLEX4_H
#define ALGPY_COMPLEX4_H

#include <stdint.h>
#include "../../vectorspaces/vec16.h"
#include "complex3.h"

#define ALGPY_COMPLEX4_N ALGPY_VEC16_N

typedef vec16 complex4;

static inline complex4
complex4_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10, double x11, double x12, double x13, double x14, double x15)
{
    complex4 result = {{x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15}};
    return result;
}

/* Write a = u + v*i3 and evaluate f(u +/- i0*v) in complex3. Reconstruct
 * u' = (p+q)/2 and v' = (p-q)/(2*i0), as in complex3_apply.
 * Principal branches are selected independently in each complex component. */
static inline complex4
complex4_apply(complex4 a, complex3 (*f)(complex3))
{
    complex3 p = f(complex3_make(
        a.data[0] - a.data[9],
        a.data[1] + a.data[8],
        a.data[2] - a.data[11],
        a.data[3] + a.data[10],
        a.data[4] - a.data[13],
        a.data[5] + a.data[12],
        a.data[6] - a.data[15],
        a.data[7] + a.data[14]));

    complex3 q = f(complex3_make(
        a.data[0] + a.data[9],
        a.data[1] - a.data[8],
        a.data[2] + a.data[11],
        a.data[3] - a.data[10],
        a.data[4] + a.data[13],
        a.data[5] - a.data[12],
        a.data[6] + a.data[15],
        a.data[7] - a.data[14]));

    return complex4_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[2] + q.data[2]) / 2,
        (p.data[3] + q.data[3]) / 2,
        (p.data[4] + q.data[4]) / 2,
        (p.data[5] + q.data[5]) / 2,
        (p.data[6] + q.data[6]) / 2,
        (p.data[7] + q.data[7]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (q.data[0] - p.data[0]) / 2,
        (p.data[3] - q.data[3]) / 2,
        (q.data[2] - p.data[2]) / 2,
        (p.data[5] - q.data[5]) / 2,
        (q.data[4] - p.data[4]) / 2,
        (p.data[7] - q.data[7]) / 2,
        (q.data[6] - p.data[6]) / 2);
}

static inline complex4
complex4_from_re(double a)
{
    return complex4_make(a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
complex4_to_re(complex4 a)
{
    return a.data[0];
}

static inline complex4
complex4_zero(void)
{
    return vec16_zero();
}

static inline complex4
complex4_one(void)
{
    return complex4_from_re(1);
}

static inline int
complex4_is_real(complex4 a)
{
    return vec16_is_real(a);
}

static inline complex4
complex4_pos(complex4 a)
{
    return vec16_pos(a);
}

static inline complex4
complex4_neg(complex4 a)
{
    return vec16_neg(a);
}

static inline complex4
complex4_re(complex4 a)
{
    return complex4_from_re(a.data[0]);
}

static inline complex4
complex4_im(complex4 a)
{
    return vec16_sub(a, complex4_re(a));
}

static inline complex4
complex4_conjugate(complex4 a)
{
    return complex4_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], -a.data[8], -a.data[9], -a.data[10], -a.data[11], -a.data[12], -a.data[13], -a.data[14], -a.data[15]);
}

static inline complex4
complex4_add(complex4 a, complex4 b)
{
    return vec16_add(a, b);
}

static inline complex4
complex4_sub(complex4 a, complex4 b)
{
    return vec16_sub(a, b);
}

static inline int
complex4_equal(complex4 a, complex4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
complex4_not_equal(complex4 a, complex4 b)
{
    return !complex4_equal(a, b);
}

static inline complex4
complex4_mul(complex4 a, complex4 b)
{
    if (complex4_is_real(b)) {
        return vec16_scale(a, b.data[0]);
    }
    if (complex4_is_real(a)) {
        return vec16_scale(b, a.data[0]);
    }
    complex3 u = complex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    complex3 v = complex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    complex3 s = complex3_make(b.data[0], b.data[1], b.data[2], b.data[3], b.data[4], b.data[5], b.data[6], b.data[7]);
    complex3 t = complex3_make(b.data[8], b.data[9], b.data[10], b.data[11], b.data[12], b.data[13], b.data[14], b.data[15]);
    complex3 lo = complex3_sub(complex3_mul(u, s), complex3_mul(v, t));
    complex3 hi = complex3_add(complex3_mul(u, t), complex3_mul(v, s));
    return complex4_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], lo.data[4], lo.data[5], lo.data[6], lo.data[7], hi.data[0], hi.data[1], hi.data[2], hi.data[3], hi.data[4], hi.data[5], hi.data[6], hi.data[7]);
}

/* Recursive norm: (u + v*i)*(u - v*i) = u*u + v*v. */
static inline complex3
complex4_norm2(complex4 a)
{
    complex3 u = complex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    complex3 v = complex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    return complex3_add(complex3_mul(u, u), complex3_mul(v, v));
}

static inline complex2
complex4_norm4(complex4 a)
{
    return complex3_norm2(complex4_norm2(a));
}

static inline complex1
complex4_norm8(complex4 a)
{
    return complex2_norm2(complex4_norm4(a));
}

static inline double
complex4_norm16(complex4 a)
{
    return complex1_norm2(complex4_norm8(a));
}

static inline double
complex4_norm(complex4 a)
{
    return sqrt(sqrt(sqrt(sqrt(complex4_norm16(a)))));
}

static inline complex4
complex4_inv(complex4 a)
{
    complex3 ninv = complex3_inv(complex4_norm2(a));
    complex3 u = complex3_make(a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7]);
    complex3 v = complex3_make(a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15]);
    complex3 lo = complex3_mul(u, ninv);
    complex3 hi = complex3_mul(complex3_neg(v), ninv);
    return complex4_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], lo.data[4], lo.data[5], lo.data[6], lo.data[7], hi.data[0], hi.data[1], hi.data[2], hi.data[3], hi.data[4], hi.data[5], hi.data[6], hi.data[7]);
}

static inline complex4
complex4_div(complex4 a, complex4 b)
{
    if (complex4_is_real(b)) {
        return vec16_div_scalar(a, b.data[0]);
    }
    return complex4_mul(a, complex4_inv(b));
}

static inline complex4
complex4_square(complex4 a)
{
    return complex4_mul(a, a);
}

static inline complex4
complex4_exp(complex4 a)
{
    return complex4_apply(a, complex3_exp);
}

static inline complex4
complex4_log(complex4 a)
{
    return complex4_apply(a, complex3_log);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline complex4
complex4_ipow(complex4 a, int64_t b)
{
    complex4 result = complex4_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return complex4_mul(a, a);
    }
    if (b == -1) {
        return complex4_inv(a);
    }
    if (b < 0) {
        a = complex4_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = complex4_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = complex4_mul(a, a);
        }
    }
    return result;
}

static inline complex4
complex4_rpow(complex4 a, complex4 b)
{
    return complex4_exp(complex4_mul(complex4_log(a), b));
}

static inline complex4
complex4_pow(complex4 a, complex4 b)
{
    return complex4_rpow(a, b);
}

static inline complex4
complex4_lpow(complex4 a, complex4 b)
{
    return complex4_rpow(a, b);
}

static inline complex4
complex4_sqrt(complex4 a)
{
    return complex4_apply(a, complex3_sqrt);
}

static inline complex4
complex4_cos(complex4 a)
{
    return complex4_apply(a, complex3_cos);
}

static inline complex4
complex4_sin(complex4 a)
{
    return complex4_apply(a, complex3_sin);
}

static inline complex4
complex4_tan(complex4 a)
{
    return complex4_apply(a, complex3_tan);
}

static inline complex4
complex4_acos(complex4 a)
{
    return complex4_apply(a, complex3_acos);
}

static inline complex4
complex4_asin(complex4 a)
{
    return complex4_apply(a, complex3_asin);
}

static inline complex4
complex4_atan(complex4 a)
{
    return complex4_apply(a, complex3_atan);
}

static inline complex4
complex4_cosh(complex4 a)
{
    return complex4_apply(a, complex3_cosh);
}

static inline complex4
complex4_sinh(complex4 a)
{
    return complex4_apply(a, complex3_sinh);
}

static inline complex4
complex4_tanh(complex4 a)
{
    return complex4_apply(a, complex3_tanh);
}

static inline complex4
complex4_acosh(complex4 a)
{
    return complex4_apply(a, complex3_acosh);
}

static inline complex4
complex4_asinh(complex4 a)
{
    return complex4_apply(a, complex3_asinh);
}

static inline complex4
complex4_atanh(complex4 a)
{
    return complex4_apply(a, complex3_atanh);
}

/* Recursive norm: (u + v*i)*(u - v*i) = u*u + v*v. */
static inline complex4
complex4_expm1(complex4 a)
{
    return complex4_apply(a, complex3_expm1);
}

static inline complex4
complex4_log1p(complex4 a)
{
    return complex4_apply(a, complex3_log1p);
}

static inline complex4
complex4_radd(double a, complex4 b)
{
    return complex4_add(complex4_from_re(a), b);
}

static inline complex4
complex4_ladd(complex4 a, double b)
{
    return complex4_add(a, complex4_from_re(b));
}

static inline complex4
complex4_rsub(double a, complex4 b)
{
    return complex4_sub(complex4_from_re(a), b);
}

static inline complex4
complex4_lsub(complex4 a, double b)
{
    return complex4_sub(a, complex4_from_re(b));
}

static inline complex4
complex4_rmul(double a, complex4 b)
{
    return vec16_scale(b, a);
}

static inline complex4
complex4_lmul(complex4 a, double b)
{
    return vec16_scale(a, b);
}

static inline complex4
complex4_rdiv(double a, complex4 b)
{
    return complex4_div(complex4_from_re(a), b);
}

static inline complex4
complex4_ldiv(complex4 a, double b)
{
    return vec16_div_scalar(a, b);
}

static inline complex4
complex4_rpows(double a, complex4 b)
{
    return complex4_rpow(complex4_from_re(a), b);
}

static inline complex4
complex4_lpows(complex4 a, double b)
{
    return complex4_rpow(a, complex4_from_re(b));
}

static inline complex4
complex4_mix(complex4 a, complex4 b, double t)
{
    complex4 delta = complex4_mul(complex4_inv(a), b);
    return complex4_mul(a, complex4_lpows(delta, t));
}

static inline void
complex4_basis(complex4 result[16])
{
    int i, j;
    for (i = 0; i < ALGPY_COMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
complex4_coeffs(complex4 result[256])
{
    complex4 basis[ALGPY_COMPLEX4_N];
    int i, j;
    complex4_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX4_N; j++) {
            result[i*ALGPY_COMPLEX4_N+j] = complex4_mul(basis[i], basis[j]);
        }
    }
}

static inline void
complex4_minpoly(complex4 a, double result[17])
{
    complex3 p = complex3_make(
        a.data[0] - a.data[9],
        a.data[1] + a.data[8],
        a.data[2] - a.data[11],
        a.data[3] + a.data[10],
        a.data[4] - a.data[13],
        a.data[5] + a.data[12],
        a.data[6] - a.data[15],
        a.data[7] + a.data[14]);
    complex3 q = complex3_make(
        a.data[0] + a.data[9],
        a.data[1] - a.data[8],
        a.data[2] + a.data[11],
        a.data[3] - a.data[10],
        a.data[4] + a.data[13],
        a.data[5] - a.data[12],
        a.data[6] + a.data[15],
        a.data[7] - a.data[14]);
    double left[9], right[9];
    int i, j;
    complex3_minpoly(p, left);
    complex3_minpoly(q, right);
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
complex4_eivals(complex4 a, complex result[16])
{
    complex3 p = complex3_make(
        a.data[0] - a.data[9],
        a.data[1] + a.data[8],
        a.data[2] - a.data[11],
        a.data[3] + a.data[10],
        a.data[4] - a.data[13],
        a.data[5] + a.data[12],
        a.data[6] - a.data[15],
        a.data[7] + a.data[14]);
    complex3 q = complex3_make(
        a.data[0] + a.data[9],
        a.data[1] - a.data[8],
        a.data[2] + a.data[11],
        a.data[3] - a.data[10],
        a.data[4] + a.data[13],
        a.data[5] - a.data[12],
        a.data[6] + a.data[15],
        a.data[7] - a.data[14]);
    complex left[8], right[8];
    complex3_eivals(p, left);
    complex3_eivals(q, right);
    /* Keep the sign-bit ordering of the generators. */
    result[0] = left[0];
    result[1] = right[1];
    result[2] = left[2];
    result[3] = right[3];
    result[4] = left[4];
    result[5] = right[5];
    result[6] = left[6];
    result[7] = right[7];
    result[8] = right[0];
    result[9] = left[1];
    result[10] = right[2];
    result[11] = left[3];
    result[12] = right[4];
    result[13] = left[5];
    result[14] = right[6];
    result[15] = left[7];
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
complex4_lmat(complex4 a, double result[ALGPY_COMPLEX4_N*ALGPY_COMPLEX4_N])
{
    complex4 basis = {{0}};
    complex4 product;
    int row, col;
    for (col = 0; col < ALGPY_COMPLEX4_N; col++) {
        basis.data[col] = 1;
        product = complex4_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_COMPLEX4_N; row++) {
            result[row * ALGPY_COMPLEX4_N + col] = product.data[row];
        }
    }
}

static inline void
complex4_rmat(complex4 a, double result[ALGPY_COMPLEX4_N*ALGPY_COMPLEX4_N])
{
    complex4 basis = {{0}};
    complex4 product;
    int row, col;
    for (row = 0; row < ALGPY_COMPLEX4_N; row++) {
        basis.data[row] = 1;
        product = complex4_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_COMPLEX4_N; col++) {
            result[row * ALGPY_COMPLEX4_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_COMPLEX4_H */
