/*
 * Shared arithmetic of bi-complex numbers w + x*i0 + y*i1 + z*i0*i1
 * with i0*i0 = i1*i1 = -1 (so (i0*i1)^2 = +1), stored as a vec4.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_COMPLEX2_H
#define ALGPY_COMPLEX2_H

#include <stdint.h>
#include "../../fields/complex.h"

#include <math.h>

#include "complex1.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_COMPLEX2_N ALGPY_VEC4_N

typedef vec4 complex2;

static inline complex2
complex2_make(double x, double x0, double x1, double x01)
{
    complex2 r = {{x, x0, x1, x01}};
    return r;
}

static inline complex2
complex2_from_re(double a)
{
    return complex2_make(a, 0, 0, 0);
}

static inline double
complex2_to_re(complex2 a)
{
    return a.data[0];
}

static inline complex2
complex2_one(void)
{
    return complex2_from_re(1);
}

static inline complex2
complex2_mul(complex2 x, complex2 y)
{
    const double *a = x.data, *b = y.data;
    return complex2_make(
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] + a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] - a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] - a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0]);
}

static inline void
complex2_basis(complex2 result[ALGPY_COMPLEX2_N])
{
    int i, j;
    for (i = 0; i < ALGPY_COMPLEX2_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX2_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
complex2_coeffs(
    complex2 result[ALGPY_COMPLEX2_N * ALGPY_COMPLEX2_N])
{
    complex2 basis[ALGPY_COMPLEX2_N];
    int i, j;
    complex2_basis(basis);
    for (i = 0; i < ALGPY_COMPLEX2_N; i++) {
        for (j = 0; j < ALGPY_COMPLEX2_N; j++) {
            result[i * ALGPY_COMPLEX2_N + j] =
                complex2_mul(basis[i], basis[j]);
        }
    }
}

/* Conjugation of the outer i1 part: (u + v*i1) -> (u - v*i1) */
static inline complex2
complex2_conjugate(complex2 a)
{
    return complex2_make(a.data[0], a.data[1], -a.data[2], -a.data[3]);
}

/* Norm a*conj(a) = u^2 + v^2 (a complex number) with a = u + v*i1. */
static inline complex
complex2_norm2(complex2 a)
{
    complex u = complex_make(a.data[0], a.data[1]);
    complex v = complex_make(a.data[2], a.data[3]);
    return vec2_add(complex_mul(u, u), complex_mul(v, v));
}

static inline double
complex2_norm4(complex2 a)
{
    return complex_norm2(complex2_norm2(a));
}

static inline double
complex2_norm(complex2 a)
{
    return sqrt(sqrt(complex2_norm4(a)));
}

static inline complex2
complex2_inv(complex2 a)
{
    complex ninv = complex_inv(complex2_norm2(a));
    complex u = complex_make(a.data[0], a.data[1]);
    complex v = complex_make(a.data[2], a.data[3]);
    complex ru = complex_mul(u, ninv);
    complex rv = complex_mul(vec2_neg(v), ninv);
    return complex2_make(ru.data[0], ru.data[1], rv.data[0], rv.data[1]);
}

static inline complex2
complex2_div(complex2 a, complex2 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return complex2_mul(a, complex2_inv(b));
}

/*
 * Apply a complex function f to a bi-complex number.
 *
 * Write a = u + v*i1 with u = w + x*i0 and v = y + z*i0 (u, v are complex
 * in i0). The algebra is commutative, so it splits into two copies of C,
 * one for each root of t^2 = -1 that i1 can map to:
 *   phi+ : i0 -> i, i1 -> +i,   a -> p = u + i*v = (w - z) + (x + y)*i
 *   phi- : i0 -> i, i1 -> -i,   a -> q = u - i*v = (w + z) + (x - y)*i
 * Then f(a) = phi+^-1(f(p)) and phi-^-1(f(q)) glued:
 *   u' = (f(p) + f(q)) / 2,  v' = (f(p) - f(q)) / (2i),
 * which gives the components computed below.
 * When f uses a principal branch, it is selected independently in p and q.
 */
static inline complex2
complex2_apply(complex2 a, complex1 (*f)(complex1))
{
    complex1 p = f(complex1_make(a.data[0] - a.data[3], a.data[1] + a.data[2]));
    complex1 q = f(complex1_make(a.data[0] + a.data[3], a.data[1] - a.data[2]));

    return complex2_make(
        (q.data[0] + p.data[0]) / 2,
        (p.data[1] + q.data[1]) / 2,
        (p.data[1] - q.data[1]) / 2,
        (q.data[0] - p.data[0]) / 2
    );
}

static inline complex2
complex2_exp(complex2 a)
{
    return complex2_apply(a, complex1_exp);
}

static inline complex2
complex2_log(complex2 a)
{
    return complex2_apply(a, complex1_log);
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline complex2
complex2_ipow(complex2 a, int64_t b)
{
    complex2 result = complex2_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return complex2_mul(a, a);
    }
    if (b == -1) {
        return complex2_inv(a);
    }
    if (b < 0) {
        a = complex2_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = complex2_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = complex2_mul(a, a);
        }
    }
    return result;
}

static inline complex2
complex2_rpow(complex2 a, complex2 b)
{
    return complex2_exp(complex2_mul(complex2_log(a), b));
}

static inline complex2
complex2_pow(complex2 a, complex2 b)
{
    return complex2_rpow(a, b);
}

static inline complex2
complex2_lpow(complex2 a, complex2 b)
{
    return complex2_rpow(a, b);
}

static inline complex2
complex2_sqrt(complex2 a)
{
    return complex2_apply(a, complex_sqrt);
}

/* Minimal polynomial coefficients (degree 4, monic), lowest degree first. */
static inline void
complex2_minpoly(complex2 a, double result[5])
{
    double x0 = a.data[0], x1 = a.data[1], x2 = a.data[2], x3 = a.data[3];
    double x02 = x0*x0, x12 = x1*x1, x22 = x2*x2, x32 = x3*x3;
    result[4] = 1;
    result[3] = -4*x0;
    result[2] = -2*x32 + 2*x12 + 2*x22 + 6*x02;
    result[1] = 4*x0*x32 - 8*x1*x2*x3 - 4*x0*x12 - 4*x0*x22 - 4*x02*x0;
    result[0] = x32*x32 + 2*x12*x32 + x12*x12 + 2*x22*x32 - 2*x02*x32
        + 8*x0*x1*x2*x3 - 2*x12*x22 + 2*x02*x12 + x22*x22 + 2*x02*x22
        + x02*x02;
}

/* Eigenvalues (x0 -/+ x3) +/- i*(x1 +/- x2) as complex numbers. */
static inline void
complex2_eivals(complex2 a, complex result[4])
{
    double x0 = a.data[0], x1 = a.data[1], x2 = a.data[2], x3 = a.data[3];
    result[0] = complex_make(x0 - x3, -x1 - x2);
    result[1] = complex_make(x0 - x3, x1 + x2);
    result[2] = complex_make(x0 + x3, x1 - x2);
    result[3] = complex_make(x0 + x3, -x1 + x2);
}

/* Common algebra API. */
static inline complex2
complex2_zero(void)
{
    return vec4_zero();
}

static inline int
complex2_is_real(complex2 a)
{
    return vec4_is_real(a);
}

static inline complex2
complex2_pos(complex2 a)
{
    return vec4_pos(a);
}

static inline complex2
complex2_neg(complex2 a)
{
    return vec4_neg(a);
}

static inline complex2
complex2_add(complex2 a, complex2 b)
{
    return vec4_add(a, b);
}

static inline complex2
complex2_sub(complex2 a, complex2 b)
{
    return vec4_sub(a, b);
}

static inline int
complex2_equal(complex2 a, complex2 b)
{
    return vec4_all_eq(a, b);
}

static inline int
complex2_not_equal(complex2 a, complex2 b)
{
    return !complex2_equal(a, b);
}

static inline complex2
complex2_re(complex2 a)
{
    return complex2_from_re(a.data[0]);
}

static inline complex2
complex2_im(complex2 a)
{
    return complex2_sub(a, complex2_re(a));
}

static inline complex2
complex2_cos(complex2 a)
{
    return complex2_apply(a, complex1_cos);
}

static inline complex2
complex2_sin(complex2 a)
{
    return complex2_apply(a, complex1_sin);
}

static inline complex2
complex2_tan(complex2 a)
{
    return complex2_apply(a, complex1_tan);
}

static inline complex2
complex2_acos(complex2 a)
{
    return complex2_apply(a, complex1_acos);
}

static inline complex2
complex2_asin(complex2 a)
{
    return complex2_apply(a, complex1_asin);
}

static inline complex2
complex2_atan(complex2 a)
{
    return complex2_apply(a, complex1_atan);
}

static inline complex2
complex2_cosh(complex2 a)
{
    return complex2_apply(a, complex1_cosh);
}

static inline complex2
complex2_sinh(complex2 a)
{
    return complex2_apply(a, complex1_sinh);
}

static inline complex2
complex2_tanh(complex2 a)
{
    return complex2_apply(a, complex1_tanh);
}

static inline complex2
complex2_acosh(complex2 a)
{
    return complex2_apply(a, complex1_acosh);
}

static inline complex2
complex2_asinh(complex2 a)
{
    return complex2_apply(a, complex1_asinh);
}

static inline complex2
complex2_atanh(complex2 a)
{
    return complex2_apply(a, complex1_atanh);
}

static inline complex2
complex2_square(complex2 a)
{
    return complex2_mul(a, a);
}

static inline complex2
complex2_radd(double a, complex2 b)
{
    return complex2_add(complex2_from_re(a), b);
}

static inline complex2
complex2_ladd(complex2 a, double b)
{
    return complex2_add(a, complex2_from_re(b));
}

static inline complex2
complex2_rsub(double a, complex2 b)
{
    return complex2_sub(complex2_from_re(a), b);
}

static inline complex2
complex2_lsub(complex2 a, double b)
{
    return complex2_sub(a, complex2_from_re(b));
}

static inline complex2
complex2_rmul(double a, complex2 b)
{
    return complex2_mul(complex2_from_re(a), b);
}

static inline complex2
complex2_lmul(complex2 a, double b)
{
    return complex2_mul(a, complex2_from_re(b));
}

static inline complex2
complex2_rdiv(double a, complex2 b)
{
    return complex2_div(complex2_from_re(a), b);
}

static inline complex2
complex2_ldiv(complex2 a, double b)
{
    return complex2_div(a, complex2_from_re(b));
}

static inline complex2
complex2_rpows(double a, complex2 b)
{
    return complex2_rpow(complex2_from_re(a), b);
}

static inline complex2
complex2_lpows(complex2 a, double b)
{
    return complex2_rpow(a, complex2_from_re(b));
}

static inline complex2
complex2_expm1(complex2 a)
{
    return complex2_lsub(complex2_exp(a), 1);
}

static inline complex2
complex2_log1p(complex2 a)
{
    return complex2_log(complex2_ladd(a, 1));
}

static inline complex2
complex2_mix(complex2 a, complex2 b, double t)
{
    complex2 delta = complex2_mul(complex2_inv(a), b);
    return complex2_mul(a, complex2_lpows(delta, t));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
complex2_lmat(complex2 a, double result[ALGPY_COMPLEX2_N*ALGPY_COMPLEX2_N])
{
    complex2 basis = {{0}};
    complex2 product;
    int row, col;
    for (col = 0; col < ALGPY_COMPLEX2_N; col++) {
        basis.data[col] = 1;
        product = complex2_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_COMPLEX2_N; row++) {
            result[row * ALGPY_COMPLEX2_N + col] = product.data[row];
        }
    }
}

static inline void
complex2_rmat(complex2 a, double result[ALGPY_COMPLEX2_N*ALGPY_COMPLEX2_N])
{
    complex2 basis = {{0}};
    complex2 product;
    int row, col;
    for (row = 0; row < ALGPY_COMPLEX2_N; row++) {
        basis.data[row] = 1;
        product = complex2_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_COMPLEX2_N; col++) {
            result[row * ALGPY_COMPLEX2_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_COMPLEX2_H */
