#ifndef ALGPY_COMPLEX3_H
#define ALGPY_COMPLEX3_H

#include "../../vectorspaces/vec8.h"
#include "complex2.h"

#define ALGPY_COMPLEX3_N ALGPY_VEC8_N

typedef vec8 complex3;

static inline complex3
complex3_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
    complex3 result = {{x0, x1, x2, x3, x4, x5, x6, x7}};
    return result;
}

/* Write a = u + v*i2 and evaluate f(u +/- i0*v) in complex2. Reconstruct
 * u' = (p+q)/2 and v' = (p-q)/(2*i0), as in complex2_apply.
 * Principal branches are selected independently in each complex component. */
static inline complex3
complex3_apply(complex3 a, complex2 (*f)(complex2))
{
    complex2 p = f(complex2_make(
        a.data[0] - a.data[5],
        a.data[1] + a.data[4],
        a.data[2] - a.data[7],
        a.data[3] + a.data[6]));

    complex2 q = f(complex2_make(
        a.data[0] + a.data[5],
        a.data[1] - a.data[4],
        a.data[2] + a.data[7],
        a.data[3] - a.data[6]));

    /* u' = (p + q)/2, v' = (p - q)/(2*i0). */
    return complex3_make(
        (p.data[0] + q.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[2] + q.data[2]) / 2,
        (p.data[3] + q.data[3]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (q.data[0] - p.data[0]) / 2,
        (p.data[3] - q.data[3]) / 2,
        (q.data[2] - p.data[2]) / 2);
}

static inline complex3
complex3_from_re(double a)
{
    return complex3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
complex3_to_re(complex3 a)
{
    return a.data[0];
}

static inline complex3
complex3_zero(void)
{
    return vec8_zero();
}

static inline complex3
complex3_one(void)
{
    return complex3_from_re(1);
}

static inline int
complex3_is_real(complex3 a)
{
    return vec8_is_real(a);
}

static inline complex3
complex3_pos(complex3 a)
{
    return vec8_pos(a);
}

static inline complex3
complex3_neg(complex3 a)
{
    return vec8_neg(a);
}

static inline complex3
complex3_re(complex3 a)
{
    return complex3_from_re(a.data[0]);
}

static inline complex3
complex3_im(complex3 a)
{
    return vec8_sub(a, complex3_re(a));
}

static inline complex3
complex3_conjugate(complex3 a)
{
    return complex3_make(a.data[0], a.data[1], a.data[2], a.data[3], -a.data[4], -a.data[5], -a.data[6], -a.data[7]);
}

static inline complex3
complex3_add(complex3 a, complex3 b)
{
    return vec8_add(a, b);
}

static inline complex3
complex3_sub(complex3 a, complex3 b)
{
    return vec8_sub(a, b);
}

static inline int
complex3_equal(complex3 a, complex3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
complex3_not_equal(complex3 a, complex3 b)
{
    return !complex3_equal(a, b);
}

static inline complex3
complex3_mul(complex3 a, complex3 b)
{
    if (complex3_is_real(b)) {
        return vec8_scale(a, b.data[0]);
    }
    if (complex3_is_real(a)) {
        return vec8_scale(b, a.data[0]);
    }
    complex2 u = complex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    complex2 v = complex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    complex2 s = complex2_make(b.data[0], b.data[1], b.data[2], b.data[3]);
    complex2 t = complex2_make(b.data[4], b.data[5], b.data[6], b.data[7]);
    complex2 lo = complex2_sub(complex2_mul(u, s), complex2_mul(v, t));
    complex2 hi = complex2_add(complex2_mul(u, t), complex2_mul(v, s));
    return complex3_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], hi.data[0], hi.data[1], hi.data[2], hi.data[3]);
}

/* Recursive norm: (u + v*i)*(u - v*i) = u*u + v*v. */
static inline complex2
complex3_norm2(complex3 a)
{
    complex2 u = complex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    complex2 v = complex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    return complex2_add(complex2_mul(u, u), complex2_mul(v, v));
}

static inline complex1
complex3_norm4(complex3 a)
{
    return complex2_norm2(complex3_norm2(a));
}

static inline double
complex3_norm8(complex3 a)
{
    return complex1_norm2(complex3_norm4(a));
}

static inline double
complex3_norm(complex3 a)
{
    return sqrt(sqrt(sqrt(complex3_norm8(a))));
}

static inline complex3
complex3_inv(complex3 a)
{
    complex2 ninv = complex2_inv(complex3_norm2(a));
    complex2 u = complex2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    complex2 v = complex2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    complex2 lo = complex2_mul(u, ninv);
    complex2 hi = complex2_mul(complex2_neg(v), ninv);
    return complex3_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], hi.data[0], hi.data[1], hi.data[2], hi.data[3]);
}

static inline complex3
complex3_div(complex3 a, complex3 b)
{
    if (complex3_is_real(b)) {
        return vec8_div_scalar(a, b.data[0]);
    }
    return complex3_mul(a, complex3_inv(b));
}

static inline complex3
complex3_square(complex3 a)
{
    return complex3_mul(a, a);
}

static inline complex3
complex3_exp(complex3 a)
{
    return complex3_apply(a, complex2_exp);
}

static inline complex3
complex3_log(complex3 a)
{
    return complex3_apply(a, complex2_log);
}

static inline complex3
complex3_pow(complex3 a, complex3 b)
{
    return complex3_exp(complex3_mul(complex3_log(a), b));
}

static inline complex3
complex3_sqrt(complex3 a)
{
    return complex3_apply(a, complex2_sqrt);
}

static inline complex3
complex3_cos(complex3 a)
{
    return complex3_apply(a, complex2_cos);
}

static inline complex3
complex3_sin(complex3 a)
{
    return complex3_apply(a, complex2_sin);
}

static inline complex3
complex3_tan(complex3 a)
{
    return complex3_apply(a, complex2_tan);
}

static inline complex3
complex3_acos(complex3 a)
{
    return complex3_apply(a, complex2_acos);
}

static inline complex3
complex3_asin(complex3 a)
{
    return complex3_apply(a, complex2_asin);
}

static inline complex3
complex3_atan(complex3 a)
{
    return complex3_apply(a, complex2_atan);
}

static inline complex3
complex3_cosh(complex3 a)
{
    return complex3_apply(a, complex2_cosh);
}

static inline complex3
complex3_sinh(complex3 a)
{
    return complex3_apply(a, complex2_sinh);
}

static inline complex3
complex3_tanh(complex3 a)
{
    return complex3_apply(a, complex2_tanh);
}

static inline complex3
complex3_acosh(complex3 a)
{
    return complex3_apply(a, complex2_acosh);
}

static inline complex3
complex3_asinh(complex3 a)
{
    return complex3_apply(a, complex2_asinh);
}

static inline complex3
complex3_atanh(complex3 a)
{
    return complex3_apply(a, complex2_atanh);
}

/* Recursive norm: (u + v*i)*(u - v*i) = u*u + v*v. */
static inline complex3
complex3_expm1(complex3 a)
{
    return complex3_apply(a, complex2_expm1);
}

static inline complex3
complex3_log1p(complex3 a)
{
    return complex3_apply(a, complex2_log1p);
}

static inline complex3
complex3_radd(double a, complex3 b)
{
    return complex3_add(complex3_from_re(a), b);
}

static inline complex3
complex3_ladd(complex3 a, double b)
{
    return complex3_add(a, complex3_from_re(b));
}

static inline complex3
complex3_rsub(double a, complex3 b)
{
    return complex3_sub(complex3_from_re(a), b);
}

static inline complex3
complex3_lsub(complex3 a, double b)
{
    return complex3_sub(a, complex3_from_re(b));
}

static inline complex3
complex3_rmul(double a, complex3 b)
{
    return vec8_scale(b, a);
}

static inline complex3
complex3_lmul(complex3 a, double b)
{
    return vec8_scale(a, b);
}

static inline complex3
complex3_rdiv(double a, complex3 b)
{
    return complex3_div(complex3_from_re(a), b);
}

static inline complex3
complex3_ldiv(complex3 a, double b)
{
    return vec8_div_scalar(a, b);
}

static inline complex3
complex3_rpow(double a, complex3 b)
{
    return complex3_pow(complex3_from_re(a), b);
}

static inline complex3
complex3_lpow(complex3 a, double b)
{
    return complex3_pow(a, complex3_from_re(b));
}

static inline complex3
complex3_mix(complex3 a, complex3 b, double t)
{
    complex3 delta = complex3_mul(complex3_inv(a), b);
    return complex3_mul(a, complex3_lpow(delta, t));
}

static inline void
complex3_basis(complex3 result[8])
{
    int i, j;
    for (i = 0; i < ALGPY_COMPLEX3_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX3_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
complex3_coeffs(complex3 result[64])
{
    complex3 basis[ALGPY_COMPLEX3_N];
    int i, j;
    complex3_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX3_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX3_N; j++) {
            result[i*ALGPY_COMPLEX3_N+j] = complex3_mul(basis[i], basis[j]);
        }
    }
}

static inline void
complex3_minpoly(complex3 a, double result[9])
{
    complex2 p = complex2_make(
        a.data[0] - a.data[5],
        a.data[1] + a.data[4],
        a.data[2] - a.data[7],
        a.data[3] + a.data[6]);
    complex2 q = complex2_make(
        a.data[0] + a.data[5],
        a.data[1] - a.data[4],
        a.data[2] + a.data[7],
        a.data[3] - a.data[6]);
    double left[5], right[5];
    int i, j;
    complex2_minpoly(p, left);
    complex2_minpoly(q, right);
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
complex3_eivals(complex3 a, complex result[8])
{
    complex2 p = complex2_make(
        a.data[0] - a.data[5],
        a.data[1] + a.data[4],
        a.data[2] - a.data[7],
        a.data[3] + a.data[6]);
    complex2 q = complex2_make(
        a.data[0] + a.data[5],
        a.data[1] - a.data[4],
        a.data[2] + a.data[7],
        a.data[3] - a.data[6]);
    complex left[4], right[4];
    complex2_eivals(p, left);
    complex2_eivals(q, right);
    /* Keep the sign-bit ordering of the generators. */
    result[0] = left[0];
    result[1] = right[2];
    result[2] = left[3];
    result[3] = right[1];
    result[4] = right[0];
    result[5] = left[2];
    result[6] = right[3];
    result[7] = left[1];
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
complex3_lmat(complex3 a, double result[ALGPY_COMPLEX3_N*ALGPY_COMPLEX3_N])
{
    complex3 coeffs[ALGPY_COMPLEX3_N * ALGPY_COMPLEX3_N];
    int row, col, i;
    complex3_coeffs(coeffs);
    for (row = 0; row < ALGPY_COMPLEX3_N; row++) {
        for (col = 0; col < ALGPY_COMPLEX3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_COMPLEX3_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_COMPLEX3_N + col].data[row];
            }
            result[row * ALGPY_COMPLEX3_N + col] = value;
        }
    }
}

static inline void
complex3_rmat(complex3 a, double result[ALGPY_COMPLEX3_N*ALGPY_COMPLEX3_N])
{
    complex3 coeffs[ALGPY_COMPLEX3_N * ALGPY_COMPLEX3_N];
    int row, col, i;
    complex3_coeffs(coeffs);
    for (row = 0; row < ALGPY_COMPLEX3_N; row++) {
        for (col = 0; col < ALGPY_COMPLEX3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_COMPLEX3_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_COMPLEX3_N + i].data[col];
            }
            result[row * ALGPY_COMPLEX3_N + col] = value;
        }
    }
}


#endif /* ALGPY_COMPLEX3_H */
