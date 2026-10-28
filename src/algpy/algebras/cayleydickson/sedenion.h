#ifndef ALGPY_SEDENION_H
#define ALGPY_SEDENION_H

#include "../general/cayleydickson/cayleydickson4.h"

#define ALGPY_SEDENION_N 16

typedef cayleydickson4 sedenion;

static inline sedenion
sedenion_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10, double x11, double x12, double x13, double x14, double x15)
{
    return cayleydickson4_make(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15);
}

static inline sedenion
sedenion_from_re(double a)
{
    return cayleydickson4_from_re(a);
}

static inline double
sedenion_to_re(sedenion a)
{
    return cayleydickson4_to_re(a);
}

static inline sedenion
sedenion_zero(void)
{
    return cayleydickson4_zero();
}

static inline sedenion
sedenion_one(void)
{
    return cayleydickson4_one();
}

static inline int
sedenion_is_real(sedenion a)
{
    return cayleydickson4_is_real(a);
}

static inline sedenion
sedenion_pos(sedenion a)
{
    return cayleydickson4_pos(a);
}

static inline sedenion
sedenion_neg(sedenion a)
{
    return cayleydickson4_neg(a);
}

static inline sedenion
sedenion_re(sedenion a)
{
    return cayleydickson4_re(a);
}

static inline sedenion
sedenion_im(sedenion a)
{
    return cayleydickson4_im(a);
}

static inline sedenion
sedenion_conjugate(sedenion a)
{
    return cayleydickson4_conjugate(a);
}

static inline sedenion
sedenion_add(sedenion a, sedenion b)
{
    return cayleydickson4_add(a, b);
}

static inline sedenion
sedenion_sub(sedenion a, sedenion b)
{
    return cayleydickson4_sub(a, b);
}

static inline int
sedenion_equal(sedenion a, sedenion b)
{
    return cayleydickson4_equal(a, b);
}

static inline int
sedenion_not_equal(sedenion a, sedenion b)
{
    return cayleydickson4_not_equal(a, b);
}

static inline sedenion
sedenion_mul(sedenion a, sedenion b)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_mul(k, a, b);
}

static inline sedenion
sedenion_div(sedenion a, sedenion b)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_div(k, a, b);
}

static inline sedenion
sedenion_pow(sedenion a, sedenion b)
{
    /* TODO: implement pow in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_pow(k, a, b);
}

static inline sedenion
sedenion_inv(sedenion a)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_inv(k, a);
}

static inline sedenion
sedenion_square(sedenion a)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_square(k, a);
}

static inline sedenion
sedenion_exp(sedenion a)
{
    /* TODO: implement exp in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_exp(k, a);
}

static inline sedenion
sedenion_log(sedenion a)
{
    /* TODO: implement log in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_log(k, a);
}

static inline sedenion
sedenion_sqrt(sedenion a)
{
    /* TODO: implement sqrt in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_sqrt(k, a);
}

static inline sedenion
sedenion_cos(sedenion a)
{
    /* TODO: implement cos in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_cos(k, a);
}

static inline sedenion
sedenion_sin(sedenion a)
{
    /* TODO: implement sin in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_sin(k, a);
}

static inline sedenion
sedenion_tan(sedenion a)
{
    /* TODO: implement tan in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_tan(k, a);
}

static inline sedenion
sedenion_acos(sedenion a)
{
    /* TODO: implement acos in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_acos(k, a);
}

static inline sedenion
sedenion_asin(sedenion a)
{
    /* TODO: implement asin in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_asin(k, a);
}

static inline sedenion
sedenion_atan(sedenion a)
{
    /* TODO: implement atan in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_atan(k, a);
}

static inline sedenion
sedenion_cosh(sedenion a)
{
    /* TODO: implement cosh in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_cosh(k, a);
}

static inline sedenion
sedenion_sinh(sedenion a)
{
    /* TODO: implement sinh in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_sinh(k, a);
}

static inline sedenion
sedenion_tanh(sedenion a)
{
    /* TODO: implement tanh in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_tanh(k, a);
}

static inline sedenion
sedenion_acosh(sedenion a)
{
    /* TODO: implement acosh in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_acosh(k, a);
}

static inline sedenion
sedenion_asinh(sedenion a)
{
    /* TODO: implement asinh in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_asinh(k, a);
}

static inline sedenion
sedenion_atanh(sedenion a)
{
    /* TODO: implement atanh in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_atanh(k, a);
}

static inline sedenion
sedenion_expm1(sedenion a)
{
    /* TODO: implement expm1 in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_expm1(k, a);
}

static inline sedenion
sedenion_log1p(sedenion a)
{
    /* TODO: implement log1p in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_log1p(k, a);
}

static inline sedenion
sedenion_radd(double a, sedenion b)
{
    return cayleydickson4_radd(a, b);
}

static inline sedenion
sedenion_ladd(sedenion a, double b)
{
    return cayleydickson4_ladd(a, b);
}

static inline sedenion
sedenion_rsub(double a, sedenion b)
{
    return cayleydickson4_rsub(a, b);
}

static inline sedenion
sedenion_lsub(sedenion a, double b)
{
    return cayleydickson4_lsub(a, b);
}

static inline sedenion
sedenion_rmul(double a, sedenion b)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_rmul(k, a, b);
}

static inline sedenion
sedenion_lmul(sedenion a, double b)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_lmul(k, a, b);
}

static inline sedenion
sedenion_rdiv(double a, sedenion b)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_rdiv(k, a, b);
}

static inline sedenion
sedenion_ldiv(sedenion a, double b)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_ldiv(k, a, b);
}

static inline sedenion
sedenion_rpow(double a, sedenion b)
{
    /* TODO: implement rpow in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_rpow(k, a, b);
}

static inline sedenion
sedenion_lpow(sedenion a, double b)
{
    /* TODO: implement lpow in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_lpow(k, a, b);
}

static inline sedenion
sedenion_mix(sedenion a, sedenion b, double t)
{
    /* TODO: implement mix in the shared general algebra. */
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_mix(k, a, b, t);
}

static inline double
sedenion_norm(sedenion a)
{
    vec4 k = {{-1, -1, -1, -1}};
    return sqrt(cayleydickson4_norm2(k, a));
}

static inline double
sedenion_norm2(sedenion a)
{
    vec4 k = {{-1, -1, -1, -1}};
    return cayleydickson4_norm2(k, a);
}

static inline void
sedenion_basis(sedenion result[16])
{
    cayleydickson4_basis(result);
}

static inline void
sedenion_coeffs(sedenion result[256])
{
    vec4 k = {{-1, -1, -1, -1}};
    cayleydickson4_coeffs(k, result);
}

static inline void
sedenion_minpoly(sedenion a, double result[3])
{
    vec4 k = {{-1, -1, -1, -1}};
    cayleydickson4_minpoly(k, a, result);
}

static inline void
sedenion_eivals(sedenion a, complex result[2])
{
    vec4 k = {{-1, -1, -1, -1}};
    cayleydickson4_eivals(k, a, result);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
sedenion_lmat(sedenion a, double result[ALGPY_SEDENION_N*ALGPY_SEDENION_N])
{
    sedenion coeffs[ALGPY_SEDENION_N * ALGPY_SEDENION_N];
    int row, col, i;
    sedenion_coeffs(coeffs);
    for (row = 0; row < ALGPY_SEDENION_N; row++) {
        for (col = 0; col < ALGPY_SEDENION_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SEDENION_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_SEDENION_N + col].data[row];
            }
            result[row * ALGPY_SEDENION_N + col] = value;
        }
    }
}

static inline void
sedenion_rmat(sedenion a, double result[ALGPY_SEDENION_N*ALGPY_SEDENION_N])
{
    sedenion coeffs[ALGPY_SEDENION_N * ALGPY_SEDENION_N];
    int row, col, i;
    sedenion_coeffs(coeffs);
    for (row = 0; row < ALGPY_SEDENION_N; row++) {
        for (col = 0; col < ALGPY_SEDENION_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_SEDENION_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_SEDENION_N + i].data[col];
            }
            result[row * ALGPY_SEDENION_N + col] = value;
        }
    }
}


#endif /* ALGPY_SEDENION_H */
