/*
 * Core arithmetic of the Cayley-Dickson algebra with three generators.
 * The parameters k are passed explicitly to every operation that depends on
 * them. This header has no Python or NumPy dependencies.
 */
#ifndef ALGPY_CAYLEYDICKSON3_H
#define ALGPY_CAYLEYDICKSON3_H

#define ALGPY_CAYLEYDICKSON3_N 8

#include <math.h>

#include "../../cayleydickson/complex.h"
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

#endif /* ALGPY_CAYLEYDICKSON3_H */
