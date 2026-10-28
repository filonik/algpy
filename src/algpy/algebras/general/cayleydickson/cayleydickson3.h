/*
 * Core arithmetic of the Cayley-Dickson algebra with three generators.
 * The parameters k are passed explicitly to every operation that depends on
 * them. This header has no Python or NumPy dependencies.
 */
#ifndef ALGPY_CAYLEYDICKSON3_H
#define ALGPY_CAYLEYDICKSON3_H

#define ALGPY_CAYLEYDICKSON3_N 8

#include <math.h>

#include "../../../fields/complex.h"
#include "cayleydickson2.h"
#include "../../../vectorspaces/vec3.h"
#include "../../../vectorspaces/vec8.h"

typedef vec8 cayleydickson3;

static inline cayleydickson3
cayleydickson3_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012)
{
    cayleydickson3 r = {{x, x0, x1, x01, x2, x02, x12, x012}};
    return r;
}

static inline cayleydickson3
cayleydickson3_re(cayleydickson3 a)
{
    return cayleydickson3_make(a.data[0], 0, 0, 0, 0, 0, 0, 0);
}

static inline cayleydickson3
cayleydickson3_im(cayleydickson3 a)
{
    return cayleydickson3_make(0, a.data[1], a.data[2], a.data[3],
                               a.data[4], a.data[5], a.data[6], a.data[7]);
}

static inline cayleydickson3
cayleydickson3_from_re(double a)
{
    return cayleydickson3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
cayleydickson3_to_re(cayleydickson3 a)
{
    return a.data[0];
}

static inline cayleydickson3
cayleydickson3_one(void)
{
    return cayleydickson3_from_re(1);
}

static inline cayleydickson3
cayleydickson3_mul(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    cayleydickson2 x0 = {{a.data[0], a.data[1], a.data[2], a.data[3]}};
    cayleydickson2 x1 = {{a.data[4], a.data[5], a.data[6], a.data[7]}};
    cayleydickson2 y0 = {{b.data[0], b.data[1], b.data[2], b.data[3]}};
    cayleydickson2 y1 = {{b.data[4], b.data[5], b.data[6], b.data[7]}};
    cayleydickson2 left = vec4_add(
        cayleydickson2_mul(k01, x0, y0),
        vec4_scale(cayleydickson2_mul(
                       k01, cayleydickson2_conjugate(y1), x1),
                   k.data[2]));
    cayleydickson2 right = vec4_add(
        cayleydickson2_mul(k01, y1, x0),
        cayleydickson2_mul(k01, x1, cayleydickson2_conjugate(y0)));
    return cayleydickson3_make(left.data[0], left.data[1],
                               left.data[2], left.data[3],
                               right.data[0], right.data[1],
                               right.data[2], right.data[3]);
}

static inline cayleydickson3
cayleydickson3_conjugate(cayleydickson3 a)
{
    return vec8_sub(cayleydickson3_re(a), cayleydickson3_im(a));
}

static inline double
cayleydickson3_norm2(vec3 k, cayleydickson3 a)
{
    double result = 0;
    int i, p, bits;
    for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        double factor = 1;
        bits = 0;
        for (p = 0; p < 3; p++) {
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

static inline cayleydickson3
cayleydickson3_inv(vec3 k, cayleydickson3 a)
{
    return vec8_div_scalar(cayleydickson3_conjugate(a),
                           cayleydickson3_norm2(k, a));
}

static inline cayleydickson3
cayleydickson3_div(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    return cayleydickson3_mul(k, a, cayleydickson3_inv(k, b));
}

static inline void
cayleydickson3_basis(cayleydickson3 result[ALGPY_CAYLEYDICKSON3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON3_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson3_coeffs(
    vec3 k,
    cayleydickson3 result[ALGPY_CAYLEYDICKSON3_N *
                          ALGPY_CAYLEYDICKSON3_N])
{
    cayleydickson3 basis[ALGPY_CAYLEYDICKSON3_N];
    int i, j;
    cayleydickson3_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON3_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON3_N + j] =
                cayleydickson3_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
cayleydickson3_minpoly(vec3 k, cayleydickson3 a, double result[3])
{
    result[0] = cayleydickson3_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

static inline void
cayleydickson3_eivals(vec3 k, cayleydickson3 a, complex result[2])
{
    complex s = complex_sqrt_real(-cayleydickson3_norm2(k, cayleydickson3_im(a)));
    result[0] = complex_make(a.data[0] - s.data[0], -s.data[1]);
    result[1] = complex_make(a.data[0] + s.data[0], +s.data[1]);
}

/* Common algebra API. Unfinished operations return NaN components. */
static inline cayleydickson3
cayleydickson3_zero(void)
{
    return vec8_zero();
}

static inline int
cayleydickson3_is_real(cayleydickson3 a)
{
    return vec8_is_real(a);
}

static inline cayleydickson3
cayleydickson3_pos(cayleydickson3 a)
{
    return vec8_pos(a);
}

static inline cayleydickson3
cayleydickson3_neg(cayleydickson3 a)
{
    return vec8_neg(a);
}

static inline cayleydickson3
cayleydickson3_add(cayleydickson3 a, cayleydickson3 b)
{
    return vec8_add(a, b);
}

static inline cayleydickson3
cayleydickson3_sub(cayleydickson3 a, cayleydickson3 b)
{
    return vec8_sub(a, b);
}

static inline int
cayleydickson3_equal(cayleydickson3 a, cayleydickson3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
cayleydickson3_not_equal(cayleydickson3 a, cayleydickson3 b)
{
    return !cayleydickson3_equal(a, b);
}

static inline cayleydickson3
cayleydickson3_exp(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_exp. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_log(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_log. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_sqrt(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_sqrt. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_cos(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_cos. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_sin(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_sin. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_tan(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_tan. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_acos(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_acos. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_asin(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_asin. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_atan(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_atan. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_cosh(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_cosh. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_sinh(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_sinh. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_tanh(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_tanh. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_acosh(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_acosh. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_asinh(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_asinh. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_atanh(vec3 k, cayleydickson3 a)
{
    /* TODO: implement cayleydickson3_atanh. */
    return cayleydickson3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline cayleydickson3
cayleydickson3_pow(vec3 k, cayleydickson3 a, cayleydickson3 b)
{
    return cayleydickson3_exp(k, cayleydickson3_mul(k, cayleydickson3_log(k, a), b));
}

static inline cayleydickson3
cayleydickson3_square(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_mul(k, a, a);
}

static inline cayleydickson3
cayleydickson3_radd(double a, cayleydickson3 b)
{
    return cayleydickson3_add(cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_ladd(cayleydickson3 a, double b)
{
    return cayleydickson3_add(a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rsub(double a, cayleydickson3 b)
{
    return cayleydickson3_sub(cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_lsub(cayleydickson3 a, double b)
{
    return cayleydickson3_sub(a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rmul(vec3 k, double a, cayleydickson3 b)
{
    return cayleydickson3_mul(k, cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_lmul(vec3 k, cayleydickson3 a, double b)
{
    return cayleydickson3_mul(k, a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rdiv(vec3 k, double a, cayleydickson3 b)
{
    return cayleydickson3_div(k, cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_ldiv(vec3 k, cayleydickson3 a, double b)
{
    return cayleydickson3_div(k, a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_rpow(vec3 k, double a, cayleydickson3 b)
{
    return cayleydickson3_pow(k, cayleydickson3_from_re(a), b);
}

static inline cayleydickson3
cayleydickson3_lpow(vec3 k, cayleydickson3 a, double b)
{
    return cayleydickson3_pow(k, a, cayleydickson3_from_re(b));
}

static inline cayleydickson3
cayleydickson3_expm1(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_lsub(cayleydickson3_exp(k, a), 1);
}

static inline cayleydickson3
cayleydickson3_log1p(vec3 k, cayleydickson3 a)
{
    return cayleydickson3_log(k, cayleydickson3_ladd(a, 1));
}

static inline cayleydickson3
cayleydickson3_mix(vec3 k, cayleydickson3 a, cayleydickson3 b, double t)
{
    cayleydickson3 delta = cayleydickson3_mul(k, cayleydickson3_inv(k, a), b);
    return cayleydickson3_mul(k, a, cayleydickson3_lpow(k, delta, t));
}

static inline complex
cayleydickson3_norm(vec3 k, cayleydickson3 a)
{
    return complex_root_real(cayleydickson3_norm2(k, a), 2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
cayleydickson3_lmat(vec3 k, cayleydickson3 a, double result[ALGPY_CAYLEYDICKSON3_N*ALGPY_CAYLEYDICKSON3_N])
{
    cayleydickson3 coeffs[ALGPY_CAYLEYDICKSON3_N * ALGPY_CAYLEYDICKSON3_N];
    int row, col, i;
    cayleydickson3_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CAYLEYDICKSON3_N; row++) {
        for (col = 0; col < ALGPY_CAYLEYDICKSON3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_CAYLEYDICKSON3_N + col].data[row];
            }
            result[row * ALGPY_CAYLEYDICKSON3_N + col] = value;
        }
    }
}

static inline void
cayleydickson3_rmat(vec3 k, cayleydickson3 a, double result[ALGPY_CAYLEYDICKSON3_N*ALGPY_CAYLEYDICKSON3_N])
{
    cayleydickson3 coeffs[ALGPY_CAYLEYDICKSON3_N * ALGPY_CAYLEYDICKSON3_N];
    int row, col, i;
    cayleydickson3_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CAYLEYDICKSON3_N; row++) {
        for (col = 0; col < ALGPY_CAYLEYDICKSON3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_CAYLEYDICKSON3_N + i].data[col];
            }
            result[row * ALGPY_CAYLEYDICKSON3_N + col] = value;
        }
    }
}


#endif /* ALGPY_CAYLEYDICKSON3_H */
