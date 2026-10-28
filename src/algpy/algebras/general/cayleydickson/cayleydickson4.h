/*
 * Core arithmetic of the Cayley-Dickson algebra with four generators.
 * The parameters k are passed explicitly to every operation that depends on
 * them. This header has no Python or NumPy dependencies.
 */
#ifndef ALGPY_CAYLEYDICKSON4_H
#define ALGPY_CAYLEYDICKSON4_H

#define ALGPY_CAYLEYDICKSON4_N 16

#include <math.h>

#include "../../../fields/complex.h"
#include "cayleydickson3.h"
#include "../../../vectorspaces/vec4.h"
#include "../../../vectorspaces/vec16.h"

typedef vec16 cayleydickson4;

static inline cayleydickson4
cayleydickson4_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012,
                    double x3, double x03, double x13, double x013,
                    double x23, double x023, double x123, double x0123)
{
    cayleydickson4 r = {{x, x0, x1, x01, x2, x02, x12, x012,
                         x3, x03, x13, x013, x23, x023, x123, x0123}};
    return r;
}

static inline cayleydickson4
cayleydickson4_re(cayleydickson4 a)
{
    return cayleydickson4_make(a.data[0], 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0);
}

static inline cayleydickson4
cayleydickson4_im(cayleydickson4 a)
{
    return cayleydickson4_make(0, a.data[1], a.data[2], a.data[3],
                               a.data[4], a.data[5], a.data[6], a.data[7],
                               a.data[8], a.data[9], a.data[10], a.data[11],
                               a.data[12], a.data[13], a.data[14], a.data[15]);
}

static inline cayleydickson4
cayleydickson4_from_re(double a)
{
    return cayleydickson4_make(a, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
cayleydickson4_to_re(cayleydickson4 a)
{
    return a.data[0];
}

static inline cayleydickson4
cayleydickson4_one(void)
{
    return cayleydickson4_from_re(1);
}

static inline cayleydickson4
cayleydickson4_mul(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    cayleydickson3 x0 = {{a.data[0], a.data[1], a.data[2], a.data[3],
                          a.data[4], a.data[5], a.data[6], a.data[7]}};
    cayleydickson3 x1 = {{a.data[8], a.data[9], a.data[10], a.data[11],
                          a.data[12], a.data[13], a.data[14], a.data[15]}};
    cayleydickson3 y0 = {{b.data[0], b.data[1], b.data[2], b.data[3],
                          b.data[4], b.data[5], b.data[6], b.data[7]}};
    cayleydickson3 y1 = {{b.data[8], b.data[9], b.data[10], b.data[11],
                          b.data[12], b.data[13], b.data[14], b.data[15]}};
    cayleydickson3 left = vec8_add(
        cayleydickson3_mul(k012, x0, y0),
        vec8_scale(cayleydickson3_mul(
                       k012, cayleydickson3_conjugate(y1), x1),
                   k.data[3]));
    cayleydickson3 right = vec8_add(
        cayleydickson3_mul(k012, y1, x0),
        cayleydickson3_mul(k012, x1, cayleydickson3_conjugate(y0)));
    return cayleydickson4_make(
        left.data[0], left.data[1], left.data[2], left.data[3],
        left.data[4], left.data[5], left.data[6], left.data[7],
        right.data[0], right.data[1], right.data[2], right.data[3],
        right.data[4], right.data[5], right.data[6], right.data[7]);
}

static inline cayleydickson4
cayleydickson4_conjugate(cayleydickson4 a)
{
    return vec16_sub(cayleydickson4_re(a), cayleydickson4_im(a));
}

static inline double
cayleydickson4_norm2(vec4 k, cayleydickson4 a)
{
    double result = 0;
    int i, p, bits;
    for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        double factor = 1;
        bits = 0;
        for (p = 0; p < 4; p++) {
            if (i & (1 << p)) {
                factor *= k.data[p];
                bits++;
            }
        }
        if (bits & 1) {
            factor = -factor;
        }
        result += factor * a.data[i] * a.data[i];
    }
    return result;
}

static inline cayleydickson4
cayleydickson4_inv(vec4 k, cayleydickson4 a)
{
    return vec16_div_scalar(cayleydickson4_conjugate(a),
                            cayleydickson4_norm2(k, a));
}

static inline cayleydickson4
cayleydickson4_div(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    return cayleydickson4_mul(k, a, cayleydickson4_inv(k, b));
}

static inline void
cayleydickson4_basis(cayleydickson4 result[ALGPY_CAYLEYDICKSON4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson4_coeffs(
    vec4 k,
    cayleydickson4 result[ALGPY_CAYLEYDICKSON4_N *
                          ALGPY_CAYLEYDICKSON4_N])
{
    cayleydickson4 basis[ALGPY_CAYLEYDICKSON4_N];
    int i, j;
    cayleydickson4_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON4_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON4_N + j] =
                cayleydickson4_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
cayleydickson4_minpoly(vec4 k, cayleydickson4 a, double result[3])
{
    result[0] = cayleydickson4_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

static inline void
cayleydickson4_eivals(vec4 k, cayleydickson4 a, complex result[2])
{
    complex s = complex_sqrt_real(-cayleydickson4_norm2(k, cayleydickson4_im(a)));
    result[0] = complex_make(a.data[0] - s.data[0], -s.data[1]);
    result[1] = complex_make(a.data[0] + s.data[0], +s.data[1]);
}

/* Common algebra API. Unfinished operations return NaN components. */
static inline cayleydickson4
cayleydickson4_zero(void)
{
    return vec16_zero();
}

static inline int
cayleydickson4_is_real(cayleydickson4 a)
{
    return vec16_is_real(a);
}

static inline cayleydickson4
cayleydickson4_pos(cayleydickson4 a)
{
    return vec16_pos(a);
}

static inline cayleydickson4
cayleydickson4_neg(cayleydickson4 a)
{
    return vec16_neg(a);
}

static inline cayleydickson4
cayleydickson4_add(cayleydickson4 a, cayleydickson4 b)
{
    return vec16_add(a, b);
}

static inline cayleydickson4
cayleydickson4_sub(cayleydickson4 a, cayleydickson4 b)
{
    return vec16_sub(a, b);
}

static inline int
cayleydickson4_equal(cayleydickson4 a, cayleydickson4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
cayleydickson4_not_equal(cayleydickson4 a, cayleydickson4 b)
{
    return !cayleydickson4_equal(a, b);
}

static inline cayleydickson4
cayleydickson4_exp(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_exp. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_log(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_log. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_sqrt(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_sqrt. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_cos(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_cos. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_sin(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_sin. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_tan(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_tan. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_acos(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_acos. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_asin(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_asin. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_atan(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_atan. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_cosh(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_cosh. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_sinh(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_sinh. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_tanh(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_tanh. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_acosh(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_acosh. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_asinh(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_asinh. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_atanh(vec4 k, cayleydickson4 a)
{
    /* TODO: implement cayleydickson4_atanh. */
    return cayleydickson4_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson4
cayleydickson4_pow(vec4 k, cayleydickson4 a, cayleydickson4 b)
{
    return cayleydickson4_exp(k, cayleydickson4_mul(k, cayleydickson4_log(k, a), b));
}

static inline cayleydickson4
cayleydickson4_square(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_mul(k, a, a);
}

static inline cayleydickson4
cayleydickson4_radd(double a, cayleydickson4 b)
{
    return cayleydickson4_add(cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_ladd(cayleydickson4 a, double b)
{
    return cayleydickson4_add(a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rsub(double a, cayleydickson4 b)
{
    return cayleydickson4_sub(cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_lsub(cayleydickson4 a, double b)
{
    return cayleydickson4_sub(a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rmul(vec4 k, double a, cayleydickson4 b)
{
    return cayleydickson4_mul(k, cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_lmul(vec4 k, cayleydickson4 a, double b)
{
    return cayleydickson4_mul(k, a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rdiv(vec4 k, double a, cayleydickson4 b)
{
    return cayleydickson4_div(k, cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_ldiv(vec4 k, cayleydickson4 a, double b)
{
    return cayleydickson4_div(k, a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_rpow(vec4 k, double a, cayleydickson4 b)
{
    return cayleydickson4_pow(k, cayleydickson4_from_re(a), b);
}

static inline cayleydickson4
cayleydickson4_lpow(vec4 k, cayleydickson4 a, double b)
{
    return cayleydickson4_pow(k, a, cayleydickson4_from_re(b));
}

static inline cayleydickson4
cayleydickson4_expm1(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_lsub(cayleydickson4_exp(k, a), 1);
}

static inline cayleydickson4
cayleydickson4_log1p(vec4 k, cayleydickson4 a)
{
    return cayleydickson4_log(k, cayleydickson4_ladd(a, 1));
}

static inline cayleydickson4
cayleydickson4_mix(vec4 k, cayleydickson4 a, cayleydickson4 b, double t)
{
    cayleydickson4 delta = cayleydickson4_mul(k, cayleydickson4_inv(k, a), b);
    return cayleydickson4_mul(k, a, cayleydickson4_lpow(k, delta, t));
}

static inline complex
cayleydickson4_norm(vec4 k, cayleydickson4 a)
{
    return complex_root_real(cayleydickson4_norm2(k, a), 2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
cayleydickson4_lmat(vec4 k, cayleydickson4 a, double result[ALGPY_CAYLEYDICKSON4_N*ALGPY_CAYLEYDICKSON4_N])
{
    cayleydickson4 coeffs[ALGPY_CAYLEYDICKSON4_N * ALGPY_CAYLEYDICKSON4_N];
    int row, col, i;
    cayleydickson4_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CAYLEYDICKSON4_N; row++) {
        for (col = 0; col < ALGPY_CAYLEYDICKSON4_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_CAYLEYDICKSON4_N + col].data[row];
            }
            result[row * ALGPY_CAYLEYDICKSON4_N + col] = value;
        }
    }
}

static inline void
cayleydickson4_rmat(vec4 k, cayleydickson4 a, double result[ALGPY_CAYLEYDICKSON4_N*ALGPY_CAYLEYDICKSON4_N])
{
    cayleydickson4 coeffs[ALGPY_CAYLEYDICKSON4_N * ALGPY_CAYLEYDICKSON4_N];
    int row, col, i;
    cayleydickson4_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CAYLEYDICKSON4_N; row++) {
        for (col = 0; col < ALGPY_CAYLEYDICKSON4_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_CAYLEYDICKSON4_N + i].data[col];
            }
            result[row * ALGPY_CAYLEYDICKSON4_N + col] = value;
        }
    }
}


#endif /* ALGPY_CAYLEYDICKSON4_H */
