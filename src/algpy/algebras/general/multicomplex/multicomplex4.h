/*
 * Core arithmetic of the multicomplex numbers
 * sum_S x_S * e_S over subsets S of {0, 1, 2, 3} (component index = bitmask of S)
 * with e0*e0 = k[0], e1*e1 = k[1], e2*e2 = k[2], e3*e3 = k[3].
 *
 * Pure C with no dependency on Python or NumPy. The parameters k are passed
 * explicitly to every operation that depends on them.
 */
#ifndef ALGPY_MULTICOMPLEX4_H
#define ALGPY_MULTICOMPLEX4_H

#define ALGPY_MULTICOMPLEX4_N 16

#include <math.h>

#include "multicomplex3.h"
#include "../../../vectorspaces/vec4.h"
#include "../../../vectorspaces/vec16.h"

typedef vec16 multicomplex4;

static inline multicomplex4
multicomplex4_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012,
                    double x3, double x03, double x13, double x013,
                    double x23, double x023, double x123, double x0123)
{
    multicomplex4 r = {{x, x0, x1, x01, x2, x02, x12, x012,
                         x3, x03, x13, x013, x23, x023, x123, x0123}};
    return r;
}

static inline multicomplex4
multicomplex4_from_re(double a)
{
    return multicomplex4_make(a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
multicomplex4_to_re(multicomplex4 a)
{
    return a.data[0];
}

static inline multicomplex4
multicomplex4_one(void)
{
    return multicomplex4_from_re(1);
}

static inline multicomplex4
multicomplex4_rsub(double a, multicomplex4 b) {
  return vec16_sub(multicomplex4_from_re(a), b);
}

static inline multicomplex4
multicomplex4_lsub(multicomplex4 a, double b) {
  return vec16_sub(a, multicomplex4_from_re(b));
}

/* e_i * e_j = (product of k[b] over the bits b shared by i and j) * e_(i^j). */
static inline multicomplex4
multicomplex4_mul(vec4 h, multicomplex4 a, multicomplex4 b)
{
    const double *k = h.data;
    multicomplex4 r;
    int i, j, c;
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        r.data[i] = 0;
    }
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX4_N; j++) {
            double t = a.data[i] * b.data[j];
            for (c = 0; c < 4; c++) {
                if ((i & j) & (1 << c)) {
                    t *= k[c];
                }
            }
            r.data[i ^ j] += t;
        }
    }
    return r;
}

static inline multicomplex4
multicomplex4_conjugate(multicomplex4 a)
{
    multicomplex4 r = a;
    int i;
    for (i = 8; i < ALGPY_MULTICOMPLEX4_N; i++) {
        r.data[i] = -r.data[i];
    }
    return r;
}

static inline multicomplex3
multicomplex4_lower(multicomplex4 a)
{
    return multicomplex3_make(a.data[0], a.data[1], a.data[2], a.data[3],
                              a.data[4], a.data[5], a.data[6], a.data[7]);
}

static inline multicomplex3
multicomplex4_upper(multicomplex4 a)
{
    return multicomplex3_make(a.data[8], a.data[9], a.data[10], a.data[11],
                              a.data[12], a.data[13], a.data[14], a.data[15]);
}

static inline multicomplex4
multicomplex4_join(multicomplex3 u, multicomplex3 v)
{
    return multicomplex4_make(u.data[0], u.data[1], u.data[2], u.data[3],
                              u.data[4], u.data[5], u.data[6], u.data[7],
                              v.data[0], v.data[1], v.data[2], v.data[3],
                              v.data[4], v.data[5], v.data[6], v.data[7]);
}

/* Norm a*conj(a) = u^2 - k[3]*v^2 with a = u + v*e3. */
static inline multicomplex3
multicomplex4_norm2(vec4 k, multicomplex4 a)
{
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    multicomplex3 u = multicomplex4_lower(a);
    multicomplex3 v = multicomplex4_upper(a);
    multicomplex3 uu = multicomplex3_mul(k012, u, u);
    multicomplex3 vv = multicomplex3_mul(k012, v, v);
    return vec8_sub(uu, vec8_scale(vv, k.data[3]));
}

static inline multicomplex2
multicomplex4_norm4(vec4 k, multicomplex4 a)
{
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    return multicomplex3_norm2(k012, multicomplex4_norm2(k, a));
}

static inline multicomplex1
multicomplex4_norm8(vec4 k, multicomplex4 a)
{
    vec2 k01 = {{k.data[0], k.data[1]}};
    return multicomplex2_norm2(k01, multicomplex4_norm4(k, a));
}

static inline double
multicomplex4_norm16(vec4 k, multicomplex4 a)
{
    vec1 k0 = {{k.data[0]}};
    return multicomplex1_norm2(k0, multicomplex4_norm8(k, a));
}

static inline complex
multicomplex4_norm(vec4 k, multicomplex4 a)
{
    return complex_root_real(multicomplex4_norm16(k, a), 16);
}

static inline multicomplex4
multicomplex4_inv(vec4 k, multicomplex4 a)
{
    vec3 k012 = {{k.data[0], k.data[1], k.data[2]}};
    multicomplex3 ninv = multicomplex3_inv(k012, multicomplex4_norm2(k, a));
    multicomplex3 u = multicomplex4_lower(a);
    multicomplex3 v = multicomplex4_upper(a);
    multicomplex3 ru = multicomplex3_mul(k012, u, ninv);
    multicomplex3 rv = multicomplex3_mul(k012, vec8_neg(v), ninv);
    return multicomplex4_join(ru, rv);
}

static inline multicomplex4
multicomplex4_div(vec4 k, multicomplex4 a, multicomplex4 b)
{
    if (vec16_is_real(b)) {
        return vec16_div_scalar(a, b.data[0]);
    }
    return multicomplex4_mul(k, a, multicomplex4_inv(k, b));
}

static inline void
multicomplex4_basis(multicomplex4 result[ALGPY_MULTICOMPLEX4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
multicomplex4_coeffs(
    vec4 k,
    multicomplex4 result[ALGPY_MULTICOMPLEX4_N * ALGPY_MULTICOMPLEX4_N])
{
    multicomplex4 basis[ALGPY_MULTICOMPLEX4_N];
    int i, j;
    multicomplex4_basis(basis);
    for (i = 0; i < ALGPY_MULTICOMPLEX4_N; i++) {
        for (j = 0; j < ALGPY_MULTICOMPLEX4_N; j++) {
            result[i * ALGPY_MULTICOMPLEX4_N + j] =
                multicomplex4_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
multicomplex4_minpoly(vec4 k, multicomplex4 a, double result[17])
{
    const int n = ALGPY_MULTICOMPLEX4_N;
    multicomplex4 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * multicomplex4_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            multicomplex4 adj = multicomplex4_lsub(uk, ck);
            uk = multicomplex4_mul(k, a, adj);
        }
    }
}

/*
 * Eigenvalues sum of x_S * prod_{i in S} (+/- sqrt(k_i)) over all subsets S,
 * with the sign of sqrt(k_i) given by the bits of the eigenvalue index,
 * as complex numbers.
 */
static inline void
multicomplex4_eivals(vec4 k, multicomplex4 a, complex result[16])
{
    complex s[4];
    complex n[16];
    int i, j, b;

    for (b = 0; b < 4; b++) {
        s[b] = complex_sqrt_real(k.data[b]);
    }
    for (i = 0; i < 16; i++) {
        complex t = complex_make(a.data[i], 0);
        for (b = 0; b < 4; b++) {
            if (i & (1 << b)) {
                t = complex_mul(t, s[b]);
            }
        }
        n[i] = t;
    }
    for (j = 0; j < 16; j++) {
        double re = 0, im = 0;
        for (i = 0; i < 16; i++) {
            int sign = 1;
            for (b = 0; b < 4; b++) {
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
