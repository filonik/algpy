/*
 * Arithmetic of real numbers, stored as a vec1.
 *
 * Pure C with no dependency on Python or NumPy.
 */
#ifndef ALGPY_REAL_H
#define ALGPY_REAL_H

#include <stdint.h>
#include "../../fields/complex.h"

#include "../../fields/double.h"
#include "../../vectorspaces/vec1.h"

#define ALGPY_REAL_N ALGPY_VEC1_N

typedef vec1 real;

static inline real
real_make(double x)
{
    real r = {{x}};
    return r;
}

static inline real
real_one(void)
{
    return real_make(1);
}

static inline real
real_mul(real x, real y)
{
    return real_make(x.data[0]*y.data[0]);
}

static inline double
real_norm2(real a)
{
    return a.data[0] * a.data[0];
}

static inline double
real_norm(real a)
{
    return fabs(a.data[0]);
}

static inline real
real_inv(real a)
{
    return real_make(1/a.data[0]);
}

static inline real
real_div(real a, real b)
{
    return real_make(a.data[0]/b.data[0]);
}

/* Common algebra API. */
static inline real
real_from_re(double a)
{
    return real_make(a);
}

static inline double
real_to_re(real a)
{
    return a.data[0];
}

static inline real
real_zero(void)
{
    return vec1_zero();
}

static inline int
real_is_real(real a)
{
    return vec1_is_real(a);
}

static inline real
real_pos(real a)
{
    return vec1_pos(a);
}

static inline real
real_neg(real a)
{
    return vec1_neg(a);
}

static inline real
real_add(real a, real b)
{
    return vec1_add(a, b);
}

static inline real
real_sub(real a, real b)
{
    return vec1_sub(a, b);
}

static inline int
real_equal(real a, real b)
{
    return vec1_all_eq(a, b);
}

static inline int
real_not_equal(real a, real b)
{
    return !real_equal(a, b);
}

static inline real
real_re(real a)
{
    return real_from_re(a.data[0]);
}

static inline real
real_im(real a)
{
    return real_sub(a, real_re(a));
}

static inline real
real_conjugate(real a)
{
    return a;
}

static inline real
real_exp(real a)
{
    return real_make(exp(a.data[0]));
}

static inline real
real_log(real a)
{
    return real_make(log(a.data[0]));
}

static inline real
real_sqrt(real a)
{
    return real_make(sqrt(a.data[0]));
}

static inline real
real_cos(real a)
{
    return real_make(cos(a.data[0]));
}

static inline real
real_sin(real a)
{
    return real_make(sin(a.data[0]));
}

static inline real
real_tan(real a)
{
    return real_make(tan(a.data[0]));
}

static inline real
real_acos(real a)
{
    return real_make(acos(a.data[0]));
}

static inline real
real_asin(real a)
{
    return real_make(asin(a.data[0]));
}

static inline real
real_atan(real a)
{
    return real_make(atan(a.data[0]));
}

static inline real
real_cosh(real a)
{
    return real_make(cosh(a.data[0]));
}

static inline real
real_sinh(real a)
{
    return real_make(sinh(a.data[0]));
}

static inline real
real_tanh(real a)
{
    return real_make(tanh(a.data[0]));
}

static inline real
real_acosh(real a)
{
    return real_make(acosh(a.data[0]));
}

static inline real
real_asinh(real a)
{
    return real_make(asinh(a.data[0]));
}

static inline real
real_atanh(real a)
{
    return real_make(atanh(a.data[0]));
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline real
real_ipow(real a, int64_t b)
{
    real result = real_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return real_mul(a, a);
    }
    if (b == -1) {
        return real_inv(a);
    }
    if (b < 0) {
        a = real_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = real_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = real_mul(a, a);
        }
    }
    return result;
}

static inline real
real_rpow(real a, real b)
{
    return real_make(double_pow(a.data[0], b.data[0]));
}

static inline real
real_pow(real a, real b)
{
    return real_rpow(a, b);
}

static inline real
real_lpow(real a, real b)
{
    return real_rpow(a, b);
}

static inline real
real_square(real a)
{
    return real_mul(a, a);
}

static inline real
real_radd(double a, real b)
{
    return real_add(real_from_re(a), b);
}

static inline real
real_ladd(real a, double b)
{
    return real_add(a, real_from_re(b));
}

static inline real
real_rsub(double a, real b)
{
    return real_sub(real_from_re(a), b);
}

static inline real
real_lsub(real a, double b)
{
    return real_sub(a, real_from_re(b));
}

static inline real
real_rmul(double a, real b)
{
    return real_mul(real_from_re(a), b);
}

static inline real
real_lmul(real a, double b)
{
    return real_mul(a, real_from_re(b));
}

static inline real
real_rdiv(double a, real b)
{
    return real_div(real_from_re(a), b);
}

static inline real
real_ldiv(real a, double b)
{
    return real_div(a, real_from_re(b));
}

static inline real
real_rpows(double a, real b)
{
    return real_rpow(real_from_re(a), b);
}

static inline real
real_lpows(real a, double b)
{
    return real_rpow(a, real_from_re(b));
}

static inline real
real_expm1(real a)
{
    return real_make(double_expm1(a.data[0]));
}

static inline real
real_log1p(real a)
{
    return real_make(double_log1p(a.data[0]));
}

static inline real
real_mix(real a, real b, double t)
{
    real delta = real_mul(real_inv(a), b);
    return real_mul(a, real_lpows(delta, t));
}

static inline void
real_basis(real result[1])
{
    result[0] = real_one();
}

static inline void
real_coeffs(real result[1])
{
    result[0] = real_one();
}

static inline void
real_minpoly(real a, double result[2])
{
    result[0] = -a.data[0];
    result[1] = 1;
}

static inline void
real_eivals(real a, complex result[1])
{
    result[0] = complex_from_re(a.data[0]);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
real_lmat(real a, double result[ALGPY_REAL_N*ALGPY_REAL_N])
{
    real basis = {{0}};
    real product;
    int row, col;
    for (col = 0; col < ALGPY_REAL_N; col++) {
        basis.data[col] = 1;
        product = real_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_REAL_N; row++) {
            result[row * ALGPY_REAL_N + col] = product.data[row];
        }
    }
}

static inline void
real_rmat(real a, double result[ALGPY_REAL_N*ALGPY_REAL_N])
{
    real basis = {{0}};
    real product;
    int row, col;
    for (row = 0; row < ALGPY_REAL_N; row++) {
        basis.data[row] = 1;
        product = real_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_REAL_N; col++) {
            result[row * ALGPY_REAL_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_REAL_H */
