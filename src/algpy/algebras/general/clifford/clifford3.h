#ifndef ALGPY_CLIFFORD3_H
#define ALGPY_CLIFFORD3_H

#define ALGPY_CLIFFORD3_N 8

#include <math.h>

#include "../../../fields/complex.h"
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

/* Common algebra API. Unfinished operations return NaN components. */
static inline clifford3
clifford3_zero(void)
{
    return vec8_zero();
}

static inline int
clifford3_is_real(clifford3 a)
{
    return vec8_is_real(a);
}

static inline clifford3
clifford3_pos(clifford3 a)
{
    return vec8_pos(a);
}

static inline clifford3
clifford3_neg(clifford3 a)
{
    return vec8_neg(a);
}

static inline clifford3
clifford3_add(clifford3 a, clifford3 b)
{
    return vec8_add(a, b);
}

static inline clifford3
clifford3_sub(clifford3 a, clifford3 b)
{
    return vec8_sub(a, b);
}

static inline int
clifford3_equal(clifford3 a, clifford3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
clifford3_not_equal(clifford3 a, clifford3 b)
{
    return !clifford3_equal(a, b);
}

static inline clifford3
clifford3_re(clifford3 a)
{
    return clifford3_from_re(a.data[0]);
}

static inline clifford3
clifford3_im(clifford3 a)
{
    return clifford3_sub(a, clifford3_re(a));
}

static inline clifford3
clifford3_exp(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_exp. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_log(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_log. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_sqrt(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_sqrt. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_cos(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_cos. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_sin(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_sin. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_tan(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_tan. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_acos(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_acos. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_asin(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_asin. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_atan(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_atan. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_cosh(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_cosh. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_sinh(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_sinh. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_tanh(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_tanh. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_acosh(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_acosh. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_asinh(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_asinh. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_atanh(vec3 k, clifford3 a)
{
    /* TODO: implement clifford3_atanh. */
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_pow(vec3 k, clifford3 a, clifford3 b)
{
    return clifford3_exp(k, clifford3_mul(k, clifford3_log(k, a), b));
}

static inline clifford3
clifford3_square(vec3 k, clifford3 a)
{
    return clifford3_mul(k, a, a);
}

static inline clifford3
clifford3_rmul(vec3 k, double a, clifford3 b)
{
    return clifford3_mul(k, clifford3_from_re(a), b);
}

static inline clifford3
clifford3_lmul(vec3 k, clifford3 a, double b)
{
    return clifford3_mul(k, a, clifford3_from_re(b));
}

static inline clifford3
clifford3_rdiv(vec3 k, double a, clifford3 b)
{
    return clifford3_div(k, clifford3_from_re(a), b);
}

static inline clifford3
clifford3_ldiv(vec3 k, clifford3 a, double b)
{
    return clifford3_div(k, a, clifford3_from_re(b));
}

static inline clifford3
clifford3_rpow(vec3 k, double a, clifford3 b)
{
    return clifford3_pow(k, clifford3_from_re(a), b);
}

static inline clifford3
clifford3_lpow(vec3 k, clifford3 a, double b)
{
    return clifford3_pow(k, a, clifford3_from_re(b));
}

static inline clifford3
clifford3_expm1(vec3 k, clifford3 a)
{
    return clifford3_lsub(clifford3_exp(k, a), 1);
}

static inline clifford3
clifford3_log1p(vec3 k, clifford3 a)
{
    return clifford3_log(k, clifford3_ladd(a, 1));
}

static inline clifford3
clifford3_mix(vec3 k, clifford3 a, clifford3 b, double t)
{
    clifford3 delta = clifford3_mul(k, clifford3_inv(k, a), b);
    return clifford3_mul(k, a, clifford3_lpow(k, delta, t));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
clifford3_lmat(vec3 k, clifford3 a, double result[ALGPY_CLIFFORD3_N*ALGPY_CLIFFORD3_N])
{
    clifford3 coeffs[ALGPY_CLIFFORD3_N * ALGPY_CLIFFORD3_N];
    int row, col, i;
    clifford3_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CLIFFORD3_N; row++) {
        for (col = 0; col < ALGPY_CLIFFORD3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_CLIFFORD3_N + col].data[row];
            }
            result[row * ALGPY_CLIFFORD3_N + col] = value;
        }
    }
}

static inline void
clifford3_rmat(vec3 k, clifford3 a, double result[ALGPY_CLIFFORD3_N*ALGPY_CLIFFORD3_N])
{
    clifford3 coeffs[ALGPY_CLIFFORD3_N * ALGPY_CLIFFORD3_N];
    int row, col, i;
    clifford3_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CLIFFORD3_N; row++) {
        for (col = 0; col < ALGPY_CLIFFORD3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_CLIFFORD3_N + i].data[col];
            }
            result[row * ALGPY_CLIFFORD3_N + col] = value;
        }
    }
}


#endif /* ALGPY_CLIFFORD3_H */
