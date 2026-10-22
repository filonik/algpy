/*
 * Core arithmetic of the multicomplex numbers
 * w + x*e1 + y*e2 + z*e1*e2 with e1*e1 = k[0], e2*e2 = k[1].
 *
 * Pure C with no dependency on Python or NumPy. The parameters k are passed
 * explicitly to every operation that depends on them.
 */
#ifndef ALGPY_MULTICOMPLEX2_H
#define ALGPY_MULTICOMPLEX2_H

#define ALGPY_MULTICOMPLEX2_N 4

#include <math.h>

#include "multicomplex1.h"
#include "../../../vectorspaces/vec1.h"
#include "../../../vectorspaces/vec2.h"
#include "../../../vectorspaces/vec4.h"

typedef vec4 multicomplex2;

static inline multicomplex2
multicomplex2_make(double x, double x0, double x1, double x01)
{
    multicomplex2 r = {{x, x0, x1, x01}};
    return r;
}

static inline multicomplex2
multicomplex2_from_re(double a)
{
    return multicomplex2_make(a, 0, 0, 0);
}

static inline double
multicomplex2_to_re(multicomplex2 a)
{
    return a.data[0];
}

static inline multicomplex2
multicomplex2_one(void)
{
    return multicomplex2_from_re(1);
}


static inline multicomplex2
multicomplex2_mul(vec2 k, multicomplex2 x, multicomplex2 y)
{
    const double *a = x.data, *b = y.data;
    if (vec4_is_real(y)) {
        return vec4_scale(x, b[0]);
    }
    if (vec4_is_real(x)) {
        return vec4_scale(y, a[0]);
    }
    return multicomplex2_make(
        a[0]*b[0] + k.data[0]*a[1]*b[1] + k.data[1]*a[2]*b[2] +
            k.data[0]*k.data[1]*a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + k.data[1]*a[2]*b[3] + k.data[1]*a[3]*b[2],
        a[0]*b[2] + k.data[0]*a[1]*b[3] + a[2]*b[0] + k.data[0]*a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0]);
}

/* Conjugation of the outer e2 part: (u + v*e2) -> (u - v*e2) */
static inline multicomplex2
multicomplex2_conjugate(multicomplex2 a)
{
    return multicomplex2_make(a.data[0], a.data[1], -a.data[2], -a.data[3]);
}


/* Norm a*conj(a) = u^2 - k[1]*v^2 with a = u + v*e2. */
static inline multicomplex1
multicomplex2_norm2(vec2 k, multicomplex2 a)
{
    vec1 k0 = {{k.data[0]}};
    multicomplex1 u = multicomplex1_make(a.data[0], a.data[1]);
    multicomplex1 v = multicomplex1_make(a.data[2], a.data[3]);
    multicomplex1 uu = multicomplex1_mul(k0, u, u);
    multicomplex1 vv = multicomplex1_mul(k0, v, v);
    return vec2_sub(uu, vec2_scale(vv, k.data[1]));
}

static inline double
multicomplex2_norm4(vec2 k, multicomplex2 a)
{
    vec1 k0 = {{k.data[0]}};
    return multicomplex1_norm2(k0, multicomplex2_norm2(k, a));
}

static inline complex
multicomplex2_norm(vec2 k, multicomplex2 a)
{
    return complex_root_real(multicomplex2_norm4(k, a), 4);
}


static inline multicomplex2
multicomplex2_inv(vec2 k, multicomplex2 a)
{
    vec1 k0 = {{k.data[0]}};
    multicomplex1 ninv = multicomplex1_inv(k0, multicomplex2_norm2(k, a));
    multicomplex1 u = multicomplex1_make(a.data[0], a.data[1]);
    multicomplex1 v = multicomplex1_make(a.data[2], a.data[3]);
    multicomplex1 ru = multicomplex1_mul(k0, u, ninv);
    multicomplex1 rv = multicomplex1_mul(k0, vec2_neg(v), ninv);
    return multicomplex2_make(ru.data[0], ru.data[1], rv.data[0], rv.data[1]);
}

static inline multicomplex2
multicomplex2_div(vec2 k, multicomplex2 a, multicomplex2 b)
{
    if (vec4_is_real(b)) {
        return vec4_div_scalar(a, b.data[0]);
    }
    return multicomplex2_mul(k, a, multicomplex2_inv(k, b));
}

static inline void
multicomplex2_basis(multicomplex2 result[ALGPY_MULTICOMPLEX2_N])
{
    result[0] = multicomplex2_make(1, 0, 0, 0);
    result[1] = multicomplex2_make(0, 1, 0, 0);
    result[2] = multicomplex2_make(0, 0, 1, 0);
    result[3] = multicomplex2_make(0, 0, 0, 1);
}

static inline void
multicomplex2_coeffs(
    vec2 k,
    multicomplex2 result[ALGPY_MULTICOMPLEX2_N * ALGPY_MULTICOMPLEX2_N])
{
    multicomplex2 basis[ALGPY_MULTICOMPLEX2_N];
    int i, j;
    multicomplex2_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX2_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX2_N; j++) {
            result[i * ALGPY_MULTICOMPLEX2_N + j] =
                multicomplex2_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
multicomplex2_minpoly(vec2 k, multicomplex2 a, double result[5])
{
    double k0 = k.data[0], k1 = k.data[1];
    double x0 = a.data[0], x1 = a.data[1], x2 = a.data[2], x3 = a.data[3];
    double x02 = x0*x0, x12 = x1*x1, x22 = x2*x2, x32 = x3*x3;
    result[4] = 1;
    result[3] = -4*x0;
    result[2] = -2*k0*k1*x32 - 2*k0*x12 - 2*k1*x22 + 6*x02;
    result[1] = 4*k0*k1*x0*x32 - 8*k0*k1*x1*x2*x3 + 4*k0*x0*x12 + 4*k1*x0*x22 - 4*x02*x0;
    result[0] = k0*k0*k1*k1*x32*x32 - 2*k0*k0*k1*x12*x32 + k0*k0*x12*x12
        - 2*k0*k1*k1*x22*x32 - 2*k0*k1*x02*x32 + 8*k0*k1*x0*x1*x2*x3
        - 2*k0*k1*x12*x22 - 2*k0*x02*x12 + k1*k1*x22*x22 - 2*k1*x02*x22
        + x02*x02;
}

/*
 * Eigenvalues x0 +/- sqrt(k0)*x1 +/- sqrt(k1)*x2 +/- sqrt(k0)*sqrt(k1)*x3,
 * as complex numbers.
 */
static inline void
multicomplex2_eivals(vec2 k, multicomplex2 a, complex result[4])
{
    complex s0 = complex_sqrt_real(k.data[0]);
    complex s1 = complex_sqrt_real(k.data[1]);
    complex t1 = complex_mul(complex_make(a.data[1], 0), s0);
    complex t2 = complex_mul(complex_make(a.data[2], 0), s1);
    complex t3 = complex_mul(complex_mul(complex_make(a.data[3], 0), s0), s1);
    result[0] = complex_make(a.data[0] - t1.data[0] - t2.data[0] + t3.data[0], -t1.data[1] - t2.data[1] + t3.data[1]);
    result[1] = complex_make(a.data[0] + t1.data[0] + t2.data[0] + t3.data[0], +t1.data[1] + t2.data[1] + t3.data[1]);
    result[2] = complex_make(a.data[0] + t1.data[0] - t2.data[0] - t3.data[0], +t1.data[1] - t2.data[1] - t3.data[1]);
    result[3] = complex_make(a.data[0] - t1.data[0] + t2.data[0] - t3.data[0], -t1.data[1] + t2.data[1] - t3.data[1]);
}

#endif /* ALGPY_MULTICOMPLEX2_H */
