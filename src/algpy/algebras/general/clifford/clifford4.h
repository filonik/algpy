#ifndef ALGPY_CLIFFORD4_H
#define ALGPY_CLIFFORD4_H

#define ALGPY_CLIFFORD4_N 16

#include <math.h>

#include "../../cayleydickson/complex.h"
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
    poly4_solve(p, result);
}

#endif /* ALGPY_CLIFFORD4_H */
