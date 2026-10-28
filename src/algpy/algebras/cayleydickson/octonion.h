#ifndef ALGPY_OCTONION_H
#define ALGPY_OCTONION_H

#include "../general/cayleydickson/cayleydickson3.h"

#define ALGPY_OCTONION_N 8

typedef cayleydickson3 octonion;

static inline octonion
octonion_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
    return cayleydickson3_make(x0, x1, x2, x3, x4, x5, x6, x7);
}

static inline octonion
octonion_from_re(double a)
{
    return cayleydickson3_from_re(a);
}

static inline double
octonion_to_re(octonion a)
{
    return cayleydickson3_to_re(a);
}

static inline octonion
octonion_zero(void)
{
    return cayleydickson3_zero();
}

static inline octonion
octonion_one(void)
{
    return cayleydickson3_one();
}

static inline int
octonion_is_real(octonion a)
{
    return cayleydickson3_is_real(a);
}

static inline octonion
octonion_pos(octonion a)
{
    return cayleydickson3_pos(a);
}

static inline octonion
octonion_neg(octonion a)
{
    return cayleydickson3_neg(a);
}

static inline octonion
octonion_re(octonion a)
{
    return cayleydickson3_re(a);
}

static inline octonion
octonion_im(octonion a)
{
    return cayleydickson3_im(a);
}

static inline octonion
octonion_conjugate(octonion a)
{
    return cayleydickson3_conjugate(a);
}

static inline octonion
octonion_add(octonion a, octonion b)
{
    return cayleydickson3_add(a, b);
}

static inline octonion
octonion_sub(octonion a, octonion b)
{
    return cayleydickson3_sub(a, b);
}

static inline int
octonion_equal(octonion a, octonion b)
{
    return cayleydickson3_equal(a, b);
}

static inline int
octonion_not_equal(octonion a, octonion b)
{
    return cayleydickson3_not_equal(a, b);
}

static inline octonion
octonion_mul(octonion a, octonion b)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_mul(k, a, b);
}

static inline octonion
octonion_div(octonion a, octonion b)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_div(k, a, b);
}

static inline octonion
octonion_pow(octonion a, octonion b)
{
    /* TODO: implement pow in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_pow(k, a, b);
}

static inline octonion
octonion_inv(octonion a)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_inv(k, a);
}

static inline octonion
octonion_square(octonion a)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_square(k, a);
}

static inline octonion
octonion_exp(octonion a)
{
    /* TODO: implement exp in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_exp(k, a);
}

static inline octonion
octonion_log(octonion a)
{
    /* TODO: implement log in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_log(k, a);
}

static inline octonion
octonion_sqrt(octonion a)
{
    /* TODO: implement sqrt in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_sqrt(k, a);
}

static inline octonion
octonion_cos(octonion a)
{
    /* TODO: implement cos in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_cos(k, a);
}

static inline octonion
octonion_sin(octonion a)
{
    /* TODO: implement sin in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_sin(k, a);
}

static inline octonion
octonion_tan(octonion a)
{
    /* TODO: implement tan in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_tan(k, a);
}

static inline octonion
octonion_acos(octonion a)
{
    /* TODO: implement acos in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_acos(k, a);
}

static inline octonion
octonion_asin(octonion a)
{
    /* TODO: implement asin in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_asin(k, a);
}

static inline octonion
octonion_atan(octonion a)
{
    /* TODO: implement atan in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_atan(k, a);
}

static inline octonion
octonion_cosh(octonion a)
{
    /* TODO: implement cosh in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_cosh(k, a);
}

static inline octonion
octonion_sinh(octonion a)
{
    /* TODO: implement sinh in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_sinh(k, a);
}

static inline octonion
octonion_tanh(octonion a)
{
    /* TODO: implement tanh in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_tanh(k, a);
}

static inline octonion
octonion_acosh(octonion a)
{
    /* TODO: implement acosh in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_acosh(k, a);
}

static inline octonion
octonion_asinh(octonion a)
{
    /* TODO: implement asinh in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_asinh(k, a);
}

static inline octonion
octonion_atanh(octonion a)
{
    /* TODO: implement atanh in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_atanh(k, a);
}

static inline octonion
octonion_expm1(octonion a)
{
    /* TODO: implement expm1 in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_expm1(k, a);
}

static inline octonion
octonion_log1p(octonion a)
{
    /* TODO: implement log1p in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_log1p(k, a);
}

static inline octonion
octonion_radd(double a, octonion b)
{
    return cayleydickson3_radd(a, b);
}

static inline octonion
octonion_ladd(octonion a, double b)
{
    return cayleydickson3_ladd(a, b);
}

static inline octonion
octonion_rsub(double a, octonion b)
{
    return cayleydickson3_rsub(a, b);
}

static inline octonion
octonion_lsub(octonion a, double b)
{
    return cayleydickson3_lsub(a, b);
}

static inline octonion
octonion_rmul(double a, octonion b)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_rmul(k, a, b);
}

static inline octonion
octonion_lmul(octonion a, double b)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_lmul(k, a, b);
}

static inline octonion
octonion_rdiv(double a, octonion b)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_rdiv(k, a, b);
}

static inline octonion
octonion_ldiv(octonion a, double b)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_ldiv(k, a, b);
}

static inline octonion
octonion_rpow(double a, octonion b)
{
    /* TODO: implement rpow in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_rpow(k, a, b);
}

static inline octonion
octonion_lpow(octonion a, double b)
{
    /* TODO: implement lpow in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_lpow(k, a, b);
}

static inline octonion
octonion_mix(octonion a, octonion b, double t)
{
    /* TODO: implement mix in the shared general algebra. */
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_mix(k, a, b, t);
}

static inline double
octonion_norm(octonion a)
{
    vec3 k = {{-1, -1, -1}};
    return sqrt(cayleydickson3_norm2(k, a));
}

static inline double
octonion_norm2(octonion a)
{
    vec3 k = {{-1, -1, -1}};
    return cayleydickson3_norm2(k, a);
}

static inline void
octonion_basis(octonion result[8])
{
    cayleydickson3_basis(result);
}

static inline void
octonion_coeffs(octonion result[64])
{
    vec3 k = {{-1, -1, -1}};
    cayleydickson3_coeffs(k, result);
}

static inline void
octonion_minpoly(octonion a, double result[3])
{
    vec3 k = {{-1, -1, -1}};
    cayleydickson3_minpoly(k, a, result);
}

static inline void
octonion_eivals(octonion a, complex result[2])
{
    vec3 k = {{-1, -1, -1}};
    cayleydickson3_eivals(k, a, result);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
octonion_lmat(octonion a, double result[ALGPY_OCTONION_N*ALGPY_OCTONION_N])
{
    octonion coeffs[ALGPY_OCTONION_N * ALGPY_OCTONION_N];
    int row, col, i;
    octonion_coeffs(coeffs);
    for (row = 0; row < ALGPY_OCTONION_N; row++) {
        for (col = 0; col < ALGPY_OCTONION_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_OCTONION_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_OCTONION_N + col].data[row];
            }
            result[row * ALGPY_OCTONION_N + col] = value;
        }
    }
}

static inline void
octonion_rmat(octonion a, double result[ALGPY_OCTONION_N*ALGPY_OCTONION_N])
{
    octonion coeffs[ALGPY_OCTONION_N * ALGPY_OCTONION_N];
    int row, col, i;
    octonion_coeffs(coeffs);
    for (row = 0; row < ALGPY_OCTONION_N; row++) {
        for (col = 0; col < ALGPY_OCTONION_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_OCTONION_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_OCTONION_N + i].data[col];
            }
            result[row * ALGPY_OCTONION_N + col] = value;
        }
    }
}


#endif /* ALGPY_OCTONION_H */
