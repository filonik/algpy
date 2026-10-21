#ifndef ALGPY_CAYLEYDICKSON2_H
#define ALGPY_CAYLEYDICKSON2_H

#define ALGPY_CAYLEYDICKSON2_N 4

#include <math.h>

#include "../../complex/complex1.h"
#include "../../../vectorspaces/vec2.h"
#include "../../../vectorspaces/vec4.h"

typedef vec4 cayleydickson2;

static inline cayleydickson2
cayleydickson2_make(double w, double x, double y, double z)
{
    cayleydickson2 r = {{w, x, y, z}};
    return r;
}

static inline cayleydickson2
cayleydickson2_re(cayleydickson2 a)
{
    return cayleydickson2_make(a.data[0], 0, 0, 0);
}

static inline cayleydickson2
cayleydickson2_im(cayleydickson2 a)
{
    return cayleydickson2_make(0, a.data[1], a.data[2], a.data[3]);
}

static inline cayleydickson2
cayleydickson2_from_re(double a)
{
    return cayleydickson2_make(a, 0, 0, 0);
}

static inline double
cayleydickson2_to_re(cayleydickson2 a)
{
    return a.data[0];
}

static inline cayleydickson2
cayleydickson2_one(void)
{
    return cayleydickson2_from_re(1);
}

static inline cayleydickson2
cayleydickson2_mul(vec2 h, cayleydickson2 a, cayleydickson2 b)
{
    const double *k = h.data, *x = a.data, *y = b.data;
    return cayleydickson2_make(
      x[0]*y[0] + k[0]*x[1]*y[1] + k[1]*x[2]*y[2] - k[0]*k[1]*x[3]*y[3],
      x[0]*y[1] + x[1]*y[0] - k[1]*x[2]*y[3] + k[1]*x[3]*y[2],
      x[0]*y[2] + k[0]*x[1]*y[3] + x[2]*y[0] - k[0]*x[3]*y[1],
      x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0]);
}

static inline void
cayleydickson2_basis(cayleydickson2 result[ALGPY_CAYLEYDICKSON2_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON2_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson2_coeffs(
    vec2 h,
    cayleydickson2 result[ALGPY_CAYLEYDICKSON2_N * ALGPY_CAYLEYDICKSON2_N])
{
    cayleydickson2 basis[ALGPY_CAYLEYDICKSON2_N];
    int i, j;
    cayleydickson2_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON2_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON2_N + j] =
                cayleydickson2_mul(h, basis[i], basis[j]);
        }
    }
}

static inline cayleydickson2
cayleydickson2_conjugate(cayleydickson2 a)
{
    return vec4_sub(cayleydickson2_re(a), cayleydickson2_im(a));
}

/* Real norm a*conj(a) = w^2 - k0*x^2 - k1*y^2 + k0*k1*z^2 */
static inline double
cayleydickson2_norm2(vec2 k, cayleydickson2 a)
{
    const double *x = a.data;
    return x[0]*x[0] - k.data[0]*x[1]*x[1] - k.data[1]*x[2]*x[2] +
        k.data[0]*k.data[1]*x[3]*x[3];
}

static inline cayleydickson2
cayleydickson2_inv(vec2 k, cayleydickson2 a)
{
    return vec4_div_scalar(cayleydickson2_conjugate(a), cayleydickson2_norm2(k, a));
}

static inline cayleydickson2
cayleydickson2_div(vec2 k, cayleydickson2 a, cayleydickson2 b)
{
    return cayleydickson2_mul(k, a, cayleydickson2_inv(k, b));
}

/* Minimal polynomial t^2 - 2*w*t + norm2(a) (coefficients in ascending order) */
static inline void
cayleydickson2_minpoly(vec2 k, cayleydickson2 a, double result[3])
{
    result[0] = cayleydickson2_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/* Roots w +/- sqrt(-norm2(im(a))) of the minimal polynomial, as complex numbers. */
static inline void
cayleydickson2_eivals(vec2 k, cayleydickson2 a, complex1 result[2])
{
    complex1 s = complex1_sqrt_real(-cayleydickson2_norm2(k, cayleydickson2_im(a)));
    result[0] = complex1_make(a.data[0] - s.data[0], -s.data[1]);
    result[1] = complex1_make(a.data[0] + s.data[0], +s.data[1]);
}

#endif /* ALGPY_CAYLEYDICKSON2_H */
