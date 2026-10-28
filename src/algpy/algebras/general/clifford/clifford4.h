#ifndef ALGPY_CLIFFORD4_H
#define ALGPY_CLIFFORD4_H

#define ALGPY_CLIFFORD4_N 16

#include <math.h>

#include "../../../fields/complex.h"
#include "../../../polynomials/poly4.h"
#include "../../../vectorspaces/vec4.h"
#include "../../../vectorspaces/vec16.h"

typedef vec16 clifford4;

static inline int popcount(unsigned int x)
{
#if defined(__clang__) || defined(__GNUC__)
    return __builtin_popcount(x);
#else
    int count = 0;
    while (x != 0) {
        x &= x - 1;  /* Clear the lowest set bit. */
        ++count;
    }
    return count;
#endif
}

static inline clifford4
clifford4_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012,
                    double x3, double x03, double x13, double x013,
                    double x23, double x023, double x123, double x0123)
{
    clifford4 r = {{x, x0, x1, x01, x2, x02, x12, x012,
                         x3, x03, x13, x013, x23, x023, x123, x0123}};
    return r;
}

static inline clifford4
clifford4_grade0(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline clifford4
clifford4_grade1(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0);
}

static inline clifford4
clifford4_grade2(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, 0, 0, x[3], 0, x[5], x[6], 0,
                          0, x[9], x[10], 0, x[12], 0, 0, 0);
}

static inline clifford4
clifford4_grade3(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, 0, 0, 0, 0, 0, 0, x[7],
                          0, 0, 0, x[11], 0, x[13], x[14], 0);
}

static inline clifford4
clifford4_grade4(clifford4 a)
{
    const double *x = a.data;
    return clifford4_make(0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, x[15]);
}

static inline clifford4
clifford4_from_re(double a)
{
    return clifford4_make(a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
clifford4_to_re(clifford4 a)
{
    return a.data[0];
}

static inline clifford4
clifford4_one(void)
{
    return clifford4_from_re(1);
}

static inline clifford4
clifford4_radd(double a, clifford4 b) {
  return vec16_add(clifford4_from_re(a), b);
}

static inline clifford4
clifford4_ladd(clifford4 a, double b) {
  return vec16_add(a, clifford4_from_re(b));
}

static inline clifford4
clifford4_rsub(double a, clifford4 b) {
  return vec16_sub(clifford4_from_re(a), b);
}

static inline clifford4
clifford4_lsub(clifford4 a, double b) {
  return vec16_sub(a, clifford4_from_re(b));
}

static inline clifford4
clifford4_mul(vec4 k, clifford4 a, clifford4 b)
{
    clifford4 r = {{0}};
    int i, j, p, q;
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD4_N; j++) {
            double factor = 1;
            int swaps = 0;
            for (p = 0; p < 4; p++) {
                if (!(j & (1 << p))) {
                    continue;
                }
                for (q = p + 1; q < 4; q++) {
                    if (i & (1 << q)) {
                        swaps++;
                    }
                }
                if (i & (1 << p)) {
                    factor *= k.data[p];
                }
            }
            if (swaps & 1) {
                factor = -factor;
            }
            r.data[i ^ j] += factor * a.data[i] * b.data[j];
        }
    }
    return r;
}

/* Clifford conjugate: grade g is scaled by (-1)^(g(g+1)/2). */
static inline clifford4
clifford4_conjugate(clifford4 a)
{
    clifford4 r = a;
    int i;
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        int g = popcount(i);
        if (((g * (g + 1)) / 2) & 1) {
            r.data[i] = -r.data[i];
        }
    }
    return r;
}

static inline void
clifford4_basis(clifford4 result[ALGPY_CLIFFORD4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        double components[ALGPY_CLIFFORD4_N] = {0};
        components[i] = 1;
        for (j = 0; j < ALGPY_CLIFFORD4_N; j++) {
            result[i].data[j] = components[j];
        }
    }
}

static inline void
clifford4_coeffs(
    vec4 k,
    clifford4 result[ALGPY_CLIFFORD4_N * ALGPY_CLIFFORD4_N])
{
    clifford4 basis[ALGPY_CLIFFORD4_N];
    int i, j;
    clifford4_basis(basis);
    for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD4_N; j++) {
            result[i * ALGPY_CLIFFORD4_N + j] =
                clifford4_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
clifford4_minpoly(vec4 k, clifford4 a, double result[5])
{
    const int n = 4;
    clifford4 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * clifford4_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            clifford4 adj = clifford4_lsub(uk, ck);
            uk = clifford4_mul(k, a, adj);
        }
    }
}

/* Norm (degree 4), the constant coefficient of the characteristic polynomial. */
static inline double
clifford4_norm4(vec4 k, clifford4 a)
{
    double p[5];
    clifford4_minpoly(k, a, p);
    return p[0];
}

static inline complex
clifford4_norm(vec4 k, clifford4 a)
{
    return complex_root_real(clifford4_norm4(k, a), 4);
}

/* Inverse from the Cayley-Hamilton theorem: p(a) = 0 with p(0) = norm2. */
static inline clifford4
clifford4_inv(vec4 k, clifford4 a)
{
    double p[5];
    clifford4 a2, a3, r;
    clifford4_minpoly(k, a, p);
    a2 = clifford4_mul(k, a, a);
    a3 = clifford4_mul(k, a2, a);
    r = vec16_add(vec16_add(a3, vec16_scale(a2, p[3])),
                  vec16_add(vec16_scale(a, p[2]), clifford4_from_re(p[1])));
    return vec16_div_scalar(vec16_neg(r), p[0]);
}

static inline clifford4
clifford4_div(vec4 k, clifford4 a, clifford4 b)
{
    return clifford4_mul(k, a, clifford4_inv(k, b));
}

static inline void
clifford4_eivals(vec4 k, clifford4 a, complex result[4])
{
    double p[5];
    clifford4_minpoly(k, a, p);
    poly4_coeffs_to_roots(p, result);
}

/* Common algebra API. Unfinished operations return NaN components. */
static inline clifford4
clifford4_zero(void)
{
    return vec16_zero();
}

static inline int
clifford4_is_real(clifford4 a)
{
    return vec16_is_real(a);
}

static inline clifford4
clifford4_pos(clifford4 a)
{
    return vec16_pos(a);
}

static inline clifford4
clifford4_neg(clifford4 a)
{
    return vec16_neg(a);
}

static inline clifford4
clifford4_add(clifford4 a, clifford4 b)
{
    return vec16_add(a, b);
}

static inline clifford4
clifford4_sub(clifford4 a, clifford4 b)
{
    return vec16_sub(a, b);
}

static inline int
clifford4_equal(clifford4 a, clifford4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
clifford4_not_equal(clifford4 a, clifford4 b)
{
    return !clifford4_equal(a, b);
}

static inline clifford4
clifford4_re(clifford4 a)
{
    return clifford4_from_re(a.data[0]);
}

static inline clifford4
clifford4_im(clifford4 a)
{
    return clifford4_sub(a, clifford4_re(a));
}

static inline clifford4
clifford4_exp(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_exp. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_log(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_log. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_sqrt(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_sqrt. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_cos(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_cos. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_sin(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_sin. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_tan(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_tan. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_acos(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_acos. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_asin(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_asin. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_atan(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_atan. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_cosh(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_cosh. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_sinh(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_sinh. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_tanh(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_tanh. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_acosh(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_acosh. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_asinh(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_asinh. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_atanh(vec4 k, clifford4 a)
{
    /* TODO: implement clifford4_atanh. */
    return clifford4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford4
clifford4_pow(vec4 k, clifford4 a, clifford4 b)
{
    return clifford4_exp(k, clifford4_mul(k, clifford4_log(k, a), b));
}

static inline clifford4
clifford4_square(vec4 k, clifford4 a)
{
    return clifford4_mul(k, a, a);
}

static inline clifford4
clifford4_rmul(vec4 k, double a, clifford4 b)
{
    return clifford4_mul(k, clifford4_from_re(a), b);
}

static inline clifford4
clifford4_lmul(vec4 k, clifford4 a, double b)
{
    return clifford4_mul(k, a, clifford4_from_re(b));
}

static inline clifford4
clifford4_rdiv(vec4 k, double a, clifford4 b)
{
    return clifford4_div(k, clifford4_from_re(a), b);
}

static inline clifford4
clifford4_ldiv(vec4 k, clifford4 a, double b)
{
    return clifford4_div(k, a, clifford4_from_re(b));
}

static inline clifford4
clifford4_rpow(vec4 k, double a, clifford4 b)
{
    return clifford4_pow(k, clifford4_from_re(a), b);
}

static inline clifford4
clifford4_lpow(vec4 k, clifford4 a, double b)
{
    return clifford4_pow(k, a, clifford4_from_re(b));
}

static inline clifford4
clifford4_expm1(vec4 k, clifford4 a)
{
    return clifford4_lsub(clifford4_exp(k, a), 1);
}

static inline clifford4
clifford4_log1p(vec4 k, clifford4 a)
{
    return clifford4_log(k, clifford4_ladd(a, 1));
}

static inline clifford4
clifford4_mix(vec4 k, clifford4 a, clifford4 b, double t)
{
    clifford4 delta = clifford4_mul(k, clifford4_inv(k, a), b);
    return clifford4_mul(k, a, clifford4_lpow(k, delta, t));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
clifford4_lmat(vec4 k, clifford4 a, double result[ALGPY_CLIFFORD4_N*ALGPY_CLIFFORD4_N])
{
    clifford4 coeffs[ALGPY_CLIFFORD4_N * ALGPY_CLIFFORD4_N];
    int row, col, i;
    clifford4_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CLIFFORD4_N; row++) {
        for (col = 0; col < ALGPY_CLIFFORD4_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_CLIFFORD4_N + col].data[row];
            }
            result[row * ALGPY_CLIFFORD4_N + col] = value;
        }
    }
}

static inline void
clifford4_rmat(vec4 k, clifford4 a, double result[ALGPY_CLIFFORD4_N*ALGPY_CLIFFORD4_N])
{
    clifford4 coeffs[ALGPY_CLIFFORD4_N * ALGPY_CLIFFORD4_N];
    int row, col, i;
    clifford4_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CLIFFORD4_N; row++) {
        for (col = 0; col < ALGPY_CLIFFORD4_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CLIFFORD4_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_CLIFFORD4_N + i].data[col];
            }
            result[row * ALGPY_CLIFFORD4_N + col] = value;
        }
    }
}


#endif /* ALGPY_CLIFFORD4_H */
