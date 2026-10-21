/*
 * Core quaternion arithmetic.
 *
 * Pure C with no dependency on Python or NumPy, so that it can be reused and
 * tested independently of the wrapper code in `_quaternion.c`.
 */
#ifndef ALGPY_QUATERNION_H
#define ALGPY_QUATERNION_H

#include <math.h>

#include "../../vectorspaces/vec2.h"
#include "../../vectorspaces/vec4.h"

#define ALGPY_QUATERNION_N ALGPY_VEC4_N

typedef vec4 quaternion;

static inline quaternion
quaternion_make(double w, double x, double y, double z)
{
    quaternion q = {{w, x, y, z}};
    return q;
}

static inline quaternion
quaternion_one(void)
{
    return quaternion_make(1, 0, 0, 0);
}

/* Hamilton product */
static inline quaternion
quaternion_mul(quaternion a, quaternion b)
{
    if (vec4_is_real(b)) {
        return vec4_scale(a, b.data[0]);
    }
    if (vec4_is_real(a)) {
        return vec4_scale(b, a.data[0]);
    }
    return quaternion_make(
        a.data[0] * b.data[0] - a.data[1] * b.data[1] -
            a.data[2] * b.data[2] - a.data[3] * b.data[3],
        a.data[0] * b.data[1] + a.data[1] * b.data[0] +
            a.data[2] * b.data[3] - a.data[3] * b.data[2],
        a.data[0] * b.data[2] - a.data[1] * b.data[3] +
            a.data[2] * b.data[0] + a.data[3] * b.data[1],
        a.data[0] * b.data[3] + a.data[1] * b.data[2] -
            a.data[2] * b.data[1] + a.data[3] * b.data[0]);
}

static inline void
quaternion_basis(quaternion result[ALGPY_QUATERNION_N])
{
    int i, j;
    for (i = 0; i < ALGPY_QUATERNION_N; i++) {
            for (j = 0; j < ALGPY_QUATERNION_N; j++) {
                result[i].data[j] = i == j ? 1 : 0;
            }
    }
}

static inline void
quaternion_coeffs(quaternion result[ALGPY_QUATERNION_N * ALGPY_QUATERNION_N])
{
    quaternion basis[ALGPY_QUATERNION_N];
    int i, j;
    quaternion_basis(basis);
    for (i = 0; i < ALGPY_QUATERNION_N; i++) {
            for (j = 0; j < ALGPY_QUATERNION_N; j++) {
                result[i * ALGPY_QUATERNION_N + j] =
                    quaternion_mul(basis[i], basis[j]);
            }
    }
}

static inline double
quaternion_norm2(quaternion a)
{
    return a.data[0] * a.data[0] + a.data[1] * a.data[1] + a.data[2] * a.data[2] + a.data[3] * a.data[3];
}

static inline double
quaternion_norm(quaternion a)
{
    return sqrt(quaternion_norm2(a));
}

static inline double
quaternion_vector_norm2(quaternion a)
{
    return a.data[1] * a.data[1] + a.data[2] * a.data[2] + a.data[3] * a.data[3];
}

static inline double
quaternion_vector_norm(quaternion a)
{
    return sqrt(quaternion_vector_norm2(a));
}

static inline void
quaternion_minpoly(quaternion a, double result[3])
{
    result[0] = quaternion_norm2(a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/*
 * Eigenvalues x0 -/+ sqrt(-|v|^2) of the quaternion, as complex numbers
 * stored in a vec2 (data[0] = real part, data[1] = imaginary part).
 * Always a conjugate pair; real quaternions give a double eigenvalue.
 */
static inline void
quaternion_eivals(quaternion a, vec2 result[2])
{
    double vn = quaternion_vector_norm(a);
    result[0] = (vec2){{a.data[0], -vn}};
    result[1] = (vec2){{a.data[0], vn}};
}

static inline quaternion
quaternion_conjugate(quaternion a)
{
    return quaternion_make(a.data[0], -a.data[1], -a.data[2], -a.data[3]);
}

static inline quaternion
quaternion_inv(quaternion a)
{
    double n = quaternion_norm2(a);
    return quaternion_make(a.data[0] / n, -a.data[1] / n, -a.data[2] / n, -a.data[3] / n);
}

static inline quaternion
quaternion_div(quaternion a, quaternion b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return quaternion_mul(a, quaternion_inv(b));
}

static inline quaternion
quaternion_exp(quaternion a)
{
    double vn = quaternion_vector_norm(a);
    double ew = exp(a.data[0]);
    double s;
    if (vn == 0) {
        return quaternion_make(ew, 0, 0, 0);
    }
    s = ew * sin(vn) / vn;
    return quaternion_make(ew * cos(vn), s * a.data[1], s * a.data[2],
                           s * a.data[3]);
}

/* Principal logarithm; for negative reals the pure part points along i. */
static inline quaternion
quaternion_log(quaternion a)
{
    double vn = quaternion_vector_norm(a);
    double f;
    if (vn == 0) {
        if (a.data[0] < 0) {
            return quaternion_make(log(-a.data[0]),
                                   3.14159265358979323846, 0, 0);
        }
        return quaternion_make(log(a.data[0]), 0, 0, 0);
    }
    f = atan2(vn, a.data[0]) / vn;
    return quaternion_make(log(sqrt(quaternion_norm2(a))),
                           f * a.data[1], f * a.data[2], f * a.data[3]);
}

#endif /* ALGPY_QUATERNION_H */
