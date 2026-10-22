#ifndef ALGPY_CLIFFORD3_H
#define ALGPY_CLIFFORD3_H

#define ALGPY_CLIFFORD3_N 8

#include <math.h>

#include "../../cayleydickson/complex.h"
#include "../../../vectorspaces/vec3.h"
#include "../../../vectorspaces/vec8.h"

typedef vec8 clifford3;

static inline clifford3
clifford3_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012)
{
    clifford3 r = {{x, x0, x1, x01, x2, x02, x12, x012}};
    return r;
}

static inline clifford3
clifford3_grade0(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(x[0], 0, 0, 0, 0, 0, 0, 0);
}

static inline clifford3
clifford3_grade1(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(0, x[1], x[2], 0, x[4], 0, 0, 0);
}

static inline clifford3
clifford3_grade2(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(0, 0, 0, x[3], 0, x[5], x[6], 0);
}

static inline clifford3
clifford3_grade3(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(0, 0, 0, 0, 0, 0, 0, x[7]);
}

static inline clifford3
clifford3_from_re(double a)
{
    return clifford3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
clifford3_to_re(clifford3 a)
{
    return a.data[0];
}

static inline clifford3
clifford3_one(void)
{
    return clifford3_from_re(1);
}

static inline clifford3
clifford3_radd(double a, clifford3 b) {
  return vec8_add(clifford3_from_re(a), b);
}

static inline clifford3
clifford3_ladd(clifford3 a, double b) {
  return vec8_add(a, clifford3_from_re(b));
}

static inline clifford3
clifford3_rsub(double a, clifford3 b) {
  return vec8_sub(clifford3_from_re(a), b);
}

static inline clifford3
clifford3_lsub(clifford3 a, double b) {
  return vec8_sub(a, clifford3_from_re(b));
}

static inline clifford3
clifford3_mul(vec3 k, clifford3 a, clifford3 b)
{
    clifford3 r = {{0}};
    int i, j, p, q;
    for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD3_N; j++) {
            double factor = 1;
            int swaps = 0;
            for (p = 0; p < 3; p++) {
                if (!(j & (1 << p))) {
                    continue;
                }
                for (q = p + 1; q < 3; q++) {
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

/* Clifford conjugate: grade 1 and grade 2 components change sign. */
static inline clifford3
clifford3_conjugate(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(x[0], -x[1], -x[2], -x[3], -x[4], -x[5], -x[6], x[7]);
}

/* Product with the conjugate: lies in the centre span{1, e0*e1*e2}. */
static inline clifford3
clifford3_conjugate_product(vec3 k, clifford3 a)
{
    return clifford3_mul(k, a, clifford3_conjugate(a));
}

/* Norm (degree 4), the constant coefficient of the minimal polynomial. */
static inline double
clifford3_norm4(vec3 k, clifford3 a)
{
    clifford3 d = clifford3_conjugate_product(k, a);
    double s = -k.data[0] * k.data[1] * k.data[2];
    return d.data[0] * d.data[0] - s * d.data[7] * d.data[7];
}

static inline complex
clifford3_norm(vec3 k, clifford3 a)
{
    return complex_root_real(clifford3_norm4(k, a), 4);
}

static inline clifford3
clifford3_inv(vec3 k, clifford3 a)
{
    clifford3 d = clifford3_conjugate_product(k, a);
    double s = -k.data[0] * k.data[1] * k.data[2];
    double n = d.data[0] * d.data[0] - s * d.data[7] * d.data[7];
    clifford3 dinv = clifford3_make(d.data[0], 0, 0, 0, 0, 0, 0, -d.data[7]);
    return vec8_div_scalar(
        clifford3_mul(k, clifford3_conjugate(a), dinv), n);
}

static inline clifford3
clifford3_div(vec3 k, clifford3 a, clifford3 b)
{
    return clifford3_mul(k, a, clifford3_inv(k, b));
}

static inline void
clifford3_basis(clifford3 result[ALGPY_CLIFFORD3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
        double components[ALGPY_CLIFFORD3_N] = {0};
        components[i] = 1;
        for (j = 0; j < ALGPY_CLIFFORD3_N; j++) {
            result[i].data[j] = components[j];
        }
    }
}

static inline void
clifford3_coeffs(
    vec3 k,
    clifford3 result[ALGPY_CLIFFORD3_N * ALGPY_CLIFFORD3_N])
{
    clifford3 basis[ALGPY_CLIFFORD3_N];
    int i, j;
    clifford3_basis(basis);
    for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD3_N; j++) {
            result[i * ALGPY_CLIFFORD3_N + j] =
                clifford3_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
clifford3_minpoly(vec3 k, clifford3 a, double result[5])
{
    const int n = 4;
    clifford3 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * clifford3_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            clifford3 adj = clifford3_lsub(uk, ck);
            uk = clifford3_mul(k, a, adj);
        }
    }
}

static inline void
clifford3_eivals(vec3 k, clifford3 a, complex result[4])
{
    /*
     * The centre is span{1, z} with z = e0*e1*e2, z*z = s. Over the centre the
     * algebra is 2x2 matrices with characteristic polynomial
     * t^2 - tau*t + delta, where tau = 2*(x + x012*z) and delta = a*conj(a).
     */
    clifford3 d = clifford3_conjugate_product(k, a);
    double s = -k.data[0] * k.data[1] * k.data[2];
    complex w = complex_sqrt_real(s);
    int m;
    for (m = 0; m < 2; m++) {
        double sign = m == 0 ? -1 : 1;
        complex tau = complex_make(
            2 * (a.data[0] + sign * a.data[7] * w.data[0]),
            2 * sign * a.data[7] * w.data[1]);
        complex delta = complex_make(
            d.data[0] + sign * d.data[7] * w.data[0],
            sign * d.data[7] * w.data[1]);
        complex tau2 = complex_mul(tau, tau);
        complex root = complex_sqrt(complex_make(
            tau2.data[0] - 4 * delta.data[0],
            tau2.data[1] - 4 * delta.data[1]));
        result[2 * m] = complex_make((tau.data[0] - root.data[0]) / 2,
                                      (tau.data[1] - root.data[1]) / 2);
        result[2 * m + 1] = complex_make((tau.data[0] + root.data[0]) / 2,
                                          (tau.data[1] + root.data[1]) / 2);
    }
}

#endif /* ALGPY_CLIFFORD3_H */
