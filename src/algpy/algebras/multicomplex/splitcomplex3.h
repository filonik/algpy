#ifndef ALGPY_SPLITCOMPLEX3_H
#define ALGPY_SPLITCOMPLEX3_H

#include "../../vectorspaces/vec8.h"
#include "splitcomplex2.h"

#define ALGPY_SPLITCOMPLEX3_N ALGPY_VEC8_N

typedef vec8 splitcomplex3;

static inline splitcomplex3
splitcomplex3_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
    splitcomplex3 result = {{x0, x1, x2, x3, x4, x5, x6, x7}};
    return result;
}

/* For a = u + v*e with e*e = +1, apply f to u+v and u-v. */
static inline splitcomplex3
splitcomplex3_apply(splitcomplex3 a, splitcomplex2 (*f)(splitcomplex2))
{
    splitcomplex2 p = splitcomplex2_make(
        a.data[0] + a.data[4],
        a.data[1] + a.data[5],
        a.data[2] + a.data[6],
        a.data[3] + a.data[7]);
    splitcomplex2 q = splitcomplex2_make(
        a.data[0] - a.data[4],
        a.data[1] - a.data[5],
        a.data[2] - a.data[6],
        a.data[3] - a.data[7]);
    p = f(p);
    q = f(q);
    return splitcomplex3_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[2] + q.data[2]) / 2,
        (p.data[3] + q.data[3]) / 2,
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (p.data[2] - q.data[2]) / 2,
        (p.data[3] - q.data[3]) / 2);
}

static inline splitcomplex3
splitcomplex3_from_re(double a)
{
    return splitcomplex3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
splitcomplex3_to_re(splitcomplex3 a)
{
    return a.data[0];
}

static inline splitcomplex3
splitcomplex3_zero(void)
{
    return vec8_zero();
}

static inline splitcomplex3
splitcomplex3_one(void)
{
    return splitcomplex3_from_re(1);
}

static inline int
splitcomplex3_is_real(splitcomplex3 a)
{
    return vec8_is_real(a);
}

static inline splitcomplex3
splitcomplex3_pos(splitcomplex3 a)
{
    return vec8_pos(a);
}

static inline splitcomplex3
splitcomplex3_neg(splitcomplex3 a)
{
    return vec8_neg(a);
}

static inline splitcomplex3
splitcomplex3_re(splitcomplex3 a)
{
    return splitcomplex3_from_re(a.data[0]);
}

static inline splitcomplex3
splitcomplex3_im(splitcomplex3 a)
{
    return vec8_sub(a, splitcomplex3_re(a));
}

static inline splitcomplex3
splitcomplex3_conjugate(splitcomplex3 a)
{
    return splitcomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3], -a.data[4], -a.data[5], -a.data[6], -a.data[7]);
}

static inline splitcomplex3
splitcomplex3_add(splitcomplex3 a, splitcomplex3 b)
{
    return vec8_add(a, b);
}

static inline splitcomplex3
splitcomplex3_sub(splitcomplex3 a, splitcomplex3 b)
{
    return vec8_sub(a, b);
}

static inline int
splitcomplex3_equal(splitcomplex3 a, splitcomplex3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
splitcomplex3_not_equal(splitcomplex3 a, splitcomplex3 b)
{
    return !splitcomplex3_equal(a, b);
}

static inline splitcomplex3
splitcomplex3_mul(splitcomplex3 a, splitcomplex3 b)
{
    if (splitcomplex3_is_real(b)) {
        return vec8_scale(a, b.data[0]);
    }
    if (splitcomplex3_is_real(a)) {
        return vec8_scale(b, a.data[0]);
    }
    splitcomplex2 u = splitcomplex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    splitcomplex2 v = splitcomplex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    splitcomplex2 s = splitcomplex2_make(b.data[0], b.data[1], b.data[2], b.data[3]);
    splitcomplex2 t = splitcomplex2_make(b.data[4], b.data[5], b.data[6], b.data[7]);
    splitcomplex2 lo = splitcomplex2_add(splitcomplex2_mul(u, s), splitcomplex2_mul(v, t));
    splitcomplex2 hi = splitcomplex2_add(splitcomplex2_mul(u, t), splitcomplex2_mul(v, s));
    return splitcomplex3_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], hi.data[0], hi.data[1], hi.data[2], hi.data[3]);
}

/* Recursive norm: (u + v*e)*(u - v*e) = u*u - v*v. */
static inline splitcomplex2
splitcomplex3_norm2(splitcomplex3 a)
{
    splitcomplex2 u = splitcomplex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    splitcomplex2 v = splitcomplex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    return splitcomplex2_sub(splitcomplex2_mul(u, u), splitcomplex2_mul(v, v));
}

static inline splitcomplex1
splitcomplex3_norm4(splitcomplex3 a)
{
    return splitcomplex2_norm2(splitcomplex3_norm2(a));
}

static inline double
splitcomplex3_norm8(splitcomplex3 a)
{
    return splitcomplex1_norm2(splitcomplex3_norm4(a));
}

static inline complex
splitcomplex3_norm(splitcomplex3 a)
{
    return complex_root_real(splitcomplex3_norm8(a), 8);
}

static inline splitcomplex3
splitcomplex3_inv(splitcomplex3 a)
{
    splitcomplex2 ninv = splitcomplex2_inv(splitcomplex3_norm2(a));
    splitcomplex2 u = splitcomplex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    splitcomplex2 v = splitcomplex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    splitcomplex2 lo = splitcomplex2_mul(u, ninv);
    splitcomplex2 hi = splitcomplex2_mul(splitcomplex2_neg(v), ninv);
    return splitcomplex3_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], hi.data[0], hi.data[1], hi.data[2], hi.data[3]);
}

static inline splitcomplex3
splitcomplex3_div(splitcomplex3 a, splitcomplex3 b)
{
    if (splitcomplex3_is_real(b)) {
        return vec8_div_scalar(a, b.data[0]);
    }
    return splitcomplex3_mul(a, splitcomplex3_inv(b));
}

static inline splitcomplex3
splitcomplex3_square(splitcomplex3 a)
{
    return splitcomplex3_mul(a, a);
}

static inline splitcomplex3
splitcomplex3_exp(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_exp);
}

static inline splitcomplex3
splitcomplex3_log(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_log);
}

static inline splitcomplex3
splitcomplex3_pow(splitcomplex3 a, splitcomplex3 b)
{
    splitcomplex2 ap = splitcomplex2_make(
        a.data[0] + a.data[4],
        a.data[1] + a.data[5],
        a.data[2] + a.data[6],
        a.data[3] + a.data[7]);
    splitcomplex2 aq = splitcomplex2_make(
        a.data[0] - a.data[4],
        a.data[1] - a.data[5],
        a.data[2] - a.data[6],
        a.data[3] - a.data[7]);
    splitcomplex2 bp = splitcomplex2_make(
        b.data[0] + b.data[4],
        b.data[1] + b.data[5],
        b.data[2] + b.data[6],
        b.data[3] + b.data[7]);
    splitcomplex2 bq = splitcomplex2_make(
        b.data[0] - b.data[4],
        b.data[1] - b.data[5],
        b.data[2] - b.data[6],
        b.data[3] - b.data[7]);
    splitcomplex2 p = splitcomplex2_pow(ap, bp);
    splitcomplex2 q = splitcomplex2_pow(aq, bq);
    return splitcomplex3_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[2] + q.data[2]) / 2,
        (p.data[3] + q.data[3]) / 2,
        (p.data[0] - q.data[0]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (p.data[2] - q.data[2]) / 2,
        (p.data[3] - q.data[3]) / 2);
}

static inline splitcomplex3
splitcomplex3_sqrt(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_sqrt);
}

static inline splitcomplex3
splitcomplex3_cos(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_cos);
}

static inline splitcomplex3
splitcomplex3_sin(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_sin);
}

static inline splitcomplex3
splitcomplex3_tan(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_tan);
}

static inline splitcomplex3
splitcomplex3_acos(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_acos);
}

static inline splitcomplex3
splitcomplex3_asin(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_asin);
}

static inline splitcomplex3
splitcomplex3_atan(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_atan);
}

static inline splitcomplex3
splitcomplex3_cosh(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_cosh);
}

static inline splitcomplex3
splitcomplex3_sinh(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_sinh);
}

static inline splitcomplex3
splitcomplex3_tanh(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_tanh);
}

static inline splitcomplex3
splitcomplex3_acosh(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_acosh);
}

static inline splitcomplex3
splitcomplex3_asinh(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_asinh);
}

static inline splitcomplex3
splitcomplex3_atanh(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_atanh);
}

static inline splitcomplex3
splitcomplex3_expm1(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_expm1);
}

static inline splitcomplex3
splitcomplex3_log1p(splitcomplex3 a)
{
    return splitcomplex3_apply(a, splitcomplex2_log1p);
}

static inline splitcomplex3
splitcomplex3_radd(double a, splitcomplex3 b)
{
    return splitcomplex3_add(splitcomplex3_from_re(a), b);
}

static inline splitcomplex3
splitcomplex3_ladd(splitcomplex3 a, double b)
{
    return splitcomplex3_add(a, splitcomplex3_from_re(b));
}

static inline splitcomplex3
splitcomplex3_rsub(double a, splitcomplex3 b)
{
    return splitcomplex3_sub(splitcomplex3_from_re(a), b);
}

static inline splitcomplex3
splitcomplex3_lsub(splitcomplex3 a, double b)
{
    return splitcomplex3_sub(a, splitcomplex3_from_re(b));
}

static inline splitcomplex3
splitcomplex3_rmul(double a, splitcomplex3 b)
{
    return vec8_scale(b, a);
}

static inline splitcomplex3
splitcomplex3_lmul(splitcomplex3 a, double b)
{
    return vec8_scale(a, b);
}

static inline splitcomplex3
splitcomplex3_rdiv(double a, splitcomplex3 b)
{
    return splitcomplex3_div(splitcomplex3_from_re(a), b);
}

static inline splitcomplex3
splitcomplex3_ldiv(splitcomplex3 a, double b)
{
    return vec8_div_scalar(a, b);
}

static inline splitcomplex3
splitcomplex3_rpow(double a, splitcomplex3 b)
{
    return splitcomplex3_pow(splitcomplex3_from_re(a), b);
}

static inline splitcomplex3
splitcomplex3_lpow(splitcomplex3 a, double b)
{
    return splitcomplex3_pow(a, splitcomplex3_from_re(b));
}

static inline splitcomplex3
splitcomplex3_mix(splitcomplex3 a, splitcomplex3 b, double t)
{
    splitcomplex3 delta = splitcomplex3_mul(splitcomplex3_inv(a), b);
    return splitcomplex3_mul(a, splitcomplex3_lpow(delta, t));
}

static inline void
splitcomplex3_basis(splitcomplex3 result[ALGPY_SPLITCOMPLEX3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_SPLITCOMPLEX3_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX3_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
splitcomplex3_coeffs(splitcomplex3 result[ALGPY_SPLITCOMPLEX3_N*ALGPY_SPLITCOMPLEX3_N])
{
    splitcomplex3 basis[ALGPY_SPLITCOMPLEX3_N];
    int i, j;
    splitcomplex3_basis(basis);
    for (i = 0; i < ALGPY_SPLITCOMPLEX3_N; i++) {
        for (j = 0; j < ALGPY_SPLITCOMPLEX3_N; j++) {
            result[i*ALGPY_SPLITCOMPLEX3_N+j] = splitcomplex3_mul(basis[i], basis[j]);
        }
    }
}

static inline void
splitcomplex3_minpoly(splitcomplex3 a, double result[9])
{
    splitcomplex2 p = splitcomplex2_make(
        a.data[0] + a.data[4],
        a.data[1] + a.data[5],
        a.data[2] + a.data[6],
        a.data[3] + a.data[7]);
    splitcomplex2 q = splitcomplex2_make(
        a.data[0] - a.data[4],
        a.data[1] - a.data[5],
        a.data[2] - a.data[6],
        a.data[3] - a.data[7]);
    double left[5], right[5];
    int i, j;
    splitcomplex2_minpoly(p, left);
    splitcomplex2_minpoly(q, right);
    for (i = 0; i <= 8; i++) {
        result[i] = 0;
    }
    for (i = 0; i <= 4; i++) {
        for (j = 0; j <= 4; j++) {
            result[i+j] += left[i]*right[j];
        }
    }
}

static inline void
splitcomplex3_eivals(splitcomplex3 a, complex result[8])
{
    splitcomplex2 p = splitcomplex2_make(
        a.data[0] + a.data[4],
        a.data[1] + a.data[5],
        a.data[2] + a.data[6],
        a.data[3] + a.data[7]);
    splitcomplex2 q = splitcomplex2_make(
        a.data[0] - a.data[4],
        a.data[1] - a.data[5],
        a.data[2] - a.data[6],
        a.data[3] - a.data[7]);
    complex left[4], right[4];
    splitcomplex2_eivals(p, left);
    splitcomplex2_eivals(q, right);
    /* Preserve the sign-bit ordering of the generators. */
    result[0] = right[3];
    result[1] = right[1];
    result[2] = right[2];
    result[3] = right[0];
    result[4] = left[3];
    result[5] = left[1];
    result[6] = left[2];
    result[7] = left[0];
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
splitcomplex3_lmat(splitcomplex3 a, double result[ALGPY_SPLITCOMPLEX3_N*ALGPY_SPLITCOMPLEX3_N])
{
    splitcomplex3 coeffs[ALGPY_SPLITCOMPLEX3_N * ALGPY_SPLITCOMPLEX3_N];
    int row, col, i;
    splitcomplex3_coeffs(coeffs);
    for (row = 0; row < ALGPY_SPLITCOMPLEX3_N; row++) {
        for (col = 0; col < ALGPY_SPLITCOMPLEX3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SPLITCOMPLEX3_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_SPLITCOMPLEX3_N + col].data[row];
            }
            result[row * ALGPY_SPLITCOMPLEX3_N + col] = value;
        }
    }
}

static inline void
splitcomplex3_rmat(splitcomplex3 a, double result[ALGPY_SPLITCOMPLEX3_N*ALGPY_SPLITCOMPLEX3_N])
{
    splitcomplex3 coeffs[ALGPY_SPLITCOMPLEX3_N * ALGPY_SPLITCOMPLEX3_N];
    int row, col, i;
    splitcomplex3_coeffs(coeffs);
    for (row = 0; row < ALGPY_SPLITCOMPLEX3_N; row++) {
        for (col = 0; col < ALGPY_SPLITCOMPLEX3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SPLITCOMPLEX3_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_SPLITCOMPLEX3_N + i].data[col];
            }
            result[row * ALGPY_SPLITCOMPLEX3_N + col] = value;
        }
    }
}


#endif /* ALGPY_SPLITCOMPLEX3_H */
