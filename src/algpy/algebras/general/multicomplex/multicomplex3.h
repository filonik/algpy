/*
 * Core arithmetic of the multicomplex numbers
 * x + x0*e0 + x1*e1 + x01*e0*e1 + x2*e2 + x02*e0*e2 + x12*e1*e2 + x012*e0*e1*e2 with e0*e0 = k[0], e1*e1 = k[1], e2*e2 = k[2].
 *
 * Pure C with no dependency on Python or NumPy. The parameters k are passed
 * explicitly to every operation that depends on them.
 */
#ifndef ALGPY_MULTICOMPLEX3_H
#define ALGPY_MULTICOMPLEX3_H

#define ALGPY_MULTICOMPLEX3_N 8

#include <math.h>

#include "multicomplex2.h"
#include "../../../vectorspaces/vec3.h"
#include "../../../vectorspaces/vec8.h"

typedef vec8 multicomplex3;

static inline multicomplex3
multicomplex3_make(double x, double x0, double x1, double x01, double x2, double x02, double x12, double x012)
{
    multicomplex3 r = {{x, x0, x1, x01, x2, x02, x12, x012}};
    return r;
}

static inline multicomplex3
multicomplex3_from_re(double a)
{
    return multicomplex3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
multicomplex3_to_re(multicomplex3 a)
{
    return a.data[0];
}

static inline multicomplex3
multicomplex3_one(void)
{
    return multicomplex3_from_re(1);
}

static inline multicomplex3
multicomplex3_rsub(double a, multicomplex3 b) {
  return vec8_sub(multicomplex3_from_re(a), b);
}

static inline multicomplex3
multicomplex3_lsub(multicomplex3 a, double b) {
  return vec8_sub(a, multicomplex3_from_re(b));
}

static inline multicomplex3
multicomplex3_mul(vec3 h, multicomplex3 a, multicomplex3 b)
{
    const double *k = h.data, *x = a.data, *y = b.data;
    return multicomplex3_make(
      x[0]*y[0] + k[0]*x[1]*y[1] + k[1]*x[2]*y[2] + k[0]*k[1]*x[3]*y[3] + k[2]*x[4]*y[4] + k[0]*k[2]*x[5]*y[5] + k[1]*k[2]*x[6]*y[6] + k[0]*k[1]*k[2]*x[7]*y[7],
      x[0]*y[1] + k[1]*k[2]*x[6]*y[7] + k[1]*k[2]*x[7]*y[6] + k[1]*x[2]*y[3] + k[1]*x[3]*y[2] + k[2]*x[4]*y[5] + k[2]*x[5]*y[4] + x[1]*y[0],
      x[0]*y[2] + k[0]*k[2]*x[5]*y[7] + k[0]*k[2]*x[7]*y[5] + k[0]*x[1]*y[3] + k[0]*x[3]*y[1] + k[2]*x[4]*y[6] + k[2]*x[6]*y[4] + x[2]*y[0],
      x[0]*y[3] + k[2]*x[4]*y[7] + k[2]*x[5]*y[6] + k[2]*x[6]*y[5] + k[2]*x[7]*y[4] + x[1]*y[2] + x[2]*y[1] + x[3]*y[0],
      x[0]*y[4] + k[0]*k[1]*x[3]*y[7] + k[0]*k[1]*x[7]*y[3] + k[0]*x[1]*y[5] + k[0]*x[5]*y[1] + k[1]*x[2]*y[6] + k[1]*x[6]*y[2] + x[4]*y[0],
      x[0]*y[5] + k[1]*x[2]*y[7] + k[1]*x[3]*y[6] + k[1]*x[6]*y[3] + k[1]*x[7]*y[2] + x[1]*y[4] + x[4]*y[1] + x[5]*y[0],
      x[0]*y[6] + k[0]*x[1]*y[7] + k[0]*x[3]*y[5] + k[0]*x[5]*y[3] + k[0]*x[7]*y[1] + x[2]*y[4] + x[4]*y[2] + x[6]*y[0],
      x[0]*y[7] + x[1]*y[6] + x[2]*y[5] + x[3]*y[4] + x[4]*y[3] + x[5]*y[2] + x[6]*y[1] + x[7]*y[0]);
}

static inline multicomplex3
multicomplex3_conjugate(multicomplex3 a)
{
    return multicomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3],
                              -a.data[4], -a.data[5], -a.data[6], -a.data[7]);
}

/* Norm a*conj(a) = u^2 - k[2]*v^2 with a = u + v*e2. */
static inline multicomplex2
multicomplex3_norm2(vec3 k, multicomplex3 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    multicomplex2 u = multicomplex2_make(a.data[0], a.data[1],
                                         a.data[2], a.data[3]);
    multicomplex2 v = multicomplex2_make(a.data[4], a.data[5],
                                         a.data[6], a.data[7]);
    multicomplex2 uu = multicomplex2_mul(k01, u, u);
    multicomplex2 vv = multicomplex2_mul(k01, v, v);
    return vec4_sub(uu, vec4_scale(vv, k.data[2]));
}

static inline multicomplex1
multicomplex3_norm4(vec3 k, multicomplex3 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    return multicomplex2_norm2(k01, multicomplex3_norm2(k, a));
}

static inline double
multicomplex3_norm8(vec3 k, multicomplex3 a)
{
    vec1 k0 = {{k.data[0]}};
    return multicomplex1_norm2(k0, multicomplex3_norm4(k, a));
}

static inline complex
multicomplex3_norm(vec3 k, multicomplex3 a)
{
    return complex_root_real(multicomplex3_norm8(k, a), 8);
}


static inline multicomplex3
multicomplex3_inv(vec3 k, multicomplex3 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    multicomplex2 ninv = multicomplex2_inv(k01, multicomplex3_norm2(k, a));
    multicomplex2 u = multicomplex2_make(a.data[0], a.data[1],
                                         a.data[2], a.data[3]);
    multicomplex2 v = multicomplex2_make(a.data[4], a.data[5],
                                         a.data[6], a.data[7]);
    multicomplex2 ru = multicomplex2_mul(k01, u, ninv);
    multicomplex2 rv = multicomplex2_mul(k01, vec4_neg(v), ninv);
    return multicomplex3_make(ru.data[0], ru.data[1], ru.data[2], ru.data[3],
                              rv.data[0], rv.data[1], rv.data[2], rv.data[3]);
}

static inline multicomplex3
multicomplex3_div(vec3 k, multicomplex3 a, multicomplex3 b)
{
    if (vec8_is_real(b)) {
        return vec8_div_scalar(a, b.data[0]);
    }
    return multicomplex3_mul(k, a, multicomplex3_inv(k, b));
}

static inline void
multicomplex3_basis(multicomplex3 result[ALGPY_MULTICOMPLEX3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_MULTICOMPLEX3_N; i++) {
        double components[ALGPY_MULTICOMPLEX3_N] = {0};
        components[i] = 1;
        for (j = 0; j < ALGPY_MULTICOMPLEX3_N; j++) {
            result[i].data[j] = components[j];
        }
    }
}

static inline void
multicomplex3_coeffs(
    vec3 k,
    multicomplex3 result[ALGPY_MULTICOMPLEX3_N * ALGPY_MULTICOMPLEX3_N])
{
    multicomplex3 basis[ALGPY_MULTICOMPLEX3_N];
    int i, j;
    multicomplex3_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX3_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX3_N; j++) {
            result[i * ALGPY_MULTICOMPLEX3_N + j] =
                multicomplex3_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
multicomplex3_minpoly(vec3 k, multicomplex3 a, double result[9])
{
    const int n = ALGPY_MULTICOMPLEX3_N;
    multicomplex3 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * multicomplex3_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            multicomplex3 adj = multicomplex3_lsub(uk, ck);
            uk = multicomplex3_mul(k, a, adj);
        }
    }
}

/*
 * Eigenvalues sum of x_S * prod_{i in S} (+/- sqrt(k_i)) over all subsets S,
 * with the sign of sqrt(k_i) given by the bits of the eigenvalue index,
 * as complex numbers.
 */
static inline void
multicomplex3_eivals(vec3 k, multicomplex3 a, complex result[8])
{
    complex s[3];
    complex n[8];
    int i, j, b;

    for (b = 0; b < 3; b++) {
        s[b] = complex_sqrt_real(k.data[b]);
    }
    for (i = 0; i < 8; i++) {
        complex t = complex_make(a.data[i], 0);
        for (b = 0; b < 3; b++) {
            if (i & (1 << b)) {
                t = complex_mul(t, s[b]);
            }
        }
        n[i] = t;
    }
    for (j = 0; j < 8; j++) {
        double re = 0, im = 0;
        for (i = 0; i < 8; i++) {
            int sign = 1;
            for (b = 0; b < 3; b++) {
                if ((i & (1 << b)) && !(j & (1 << b))) {
                    sign = -sign;
                }
            }
            re += sign * n[i].data[0];
            im += sign * n[i].data[1];
        }
            result[j] = complex_make(re, im);
    }
}

#endif
