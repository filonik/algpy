#ifndef ALGPY_DUAL5_H
#define ALGPY_DUAL5_H

#include <string.h>
#include <stdint.h>
#include "../../vectorspaces/vec32.h"
#include "dual4.h"

#define ALGPY_DUAL5_N ALGPY_VEC32_N

typedef vec32 dual5;

static inline dual5
dual5_make(const double coeffs[ALGPY_DUAL5_N])
{
    dual5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

/* For a = u + v*e with e*e = 0, f(a) = f(u) + v*f'(u)*e.
 * Both callbacks are lower-dimensional dual functions, carrying all mixed
 * derivatives. They must be defined and differentiable at the real part. */
static inline dual5
dual5_apply(dual5 a, dual4 (*f)(dual4), dual4 (*deriv1)(dual4))
{
    dual4 u = dual4_make(a.data);
    dual4 v = dual4_make(a.data + 16);
    dual4 value = f(u);
    dual4 tangent = dual4_mul(v, deriv1(u));
    dual5 result;
    memcpy(result.data, value.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_DUAL5_N/2, tangent.data, sizeof(result.data)/2);
    return result;
}

static inline dual5
dual5_from_re(double a)
{
    dual5 result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
dual5_to_re(dual5 a)
{
    return a.data[0];
}

static inline dual5
dual5_zero(void)
{
    return vec32_zero();
}

static inline dual5
dual5_one(void)
{
    return dual5_from_re(1);
}

static inline int
dual5_is_real(dual5 a)
{
    return vec32_is_real(a);
}

static inline dual5
dual5_pos(dual5 a)
{
    return vec32_pos(a);
}

static inline dual5
dual5_neg(dual5 a)
{
    return vec32_neg(a);
}

static inline dual5
dual5_re(dual5 a)
{
    return dual5_from_re(a.data[0]);
}

static inline dual5
dual5_im(dual5 a)
{
    return vec32_sub(a, dual5_re(a));
}

static inline dual5
dual5_conjugate(dual5 a)
{
    return dual5_make((const double[]){a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], a.data[8], a.data[9], a.data[10], a.data[11], a.data[12], a.data[13], a.data[14], a.data[15], -a.data[16], -a.data[17], -a.data[18], -a.data[19], -a.data[20], -a.data[21], -a.data[22], -a.data[23], -a.data[24], -a.data[25], -a.data[26], -a.data[27], -a.data[28], -a.data[29], -a.data[30], -a.data[31]});
}

static inline dual5
dual5_add(dual5 a, dual5 b)
{
    return vec32_add(a, b);
}

static inline dual5
dual5_sub(dual5 a, dual5 b)
{
    return vec32_sub(a, b);
}

static inline int
dual5_equal(dual5 a, dual5 b)
{
    return vec32_all_eq(a, b);
}

static inline int
dual5_not_equal(dual5 a, dual5 b)
{
    return !dual5_equal(a, b);
}

static inline dual5
dual5_mul(dual5 a, dual5 b)
{
    if (dual5_is_real(b)) {
        return vec32_scale(a, b.data[0]);
    }
    if (dual5_is_real(a)) {
        return vec32_scale(b, a.data[0]);
    }
    dual4 u = dual4_make(a.data);
    dual4 v = dual4_make(a.data + 16);
    dual4 s = dual4_make(b.data);
    dual4 t = dual4_make(b.data + 16);
    dual4 lo = dual4_mul(u, s);
    dual4 hi = dual4_add(dual4_mul(u, t), dual4_mul(v, s));
    dual5 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_DUAL5_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

/* Norm a*conjugate(a) = u*u, since the outer generator squares to zero. */
static inline dual4
dual5_norm2(dual5 a)
{
    dual4 u = dual4_make(a.data);
    return dual4_square(u);
}

static inline dual3
dual5_norm4(dual5 a)
{
    return dual4_norm2(dual5_norm2(a));
}

static inline dual2
dual5_norm8(dual5 a)
{
    return dual3_norm2(dual5_norm4(a));
}

static inline dual1
dual5_norm16(dual5 a)
{
    return dual2_norm2(dual5_norm8(a));
}

static inline double
dual5_norm32(dual5 a)
{
    return dual1_norm2(dual5_norm16(a));
}

static inline double
dual5_norm(dual5 a)
{
    return fabs(a.data[0]);
}

static inline dual5
dual5_inv(dual5 a)
{
    dual4 u = dual4_make(a.data);
    dual4 v = dual4_make(a.data + 16);
    dual4 lo = dual4_inv(u);
    dual4 hi = dual4_neg(dual4_mul(dual4_mul(v, lo), lo));
    dual5 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_DUAL5_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

static inline dual5
dual5_div(dual5 a, dual5 b)
{
    if (dual5_is_real(b)) {
        return vec32_div_scalar(a, b.data[0]);
    }
    return dual5_mul(a, dual5_inv(b));
}

static inline dual5
dual5_square(dual5 a)
{
    return dual5_mul(a, a);
}

static inline dual5
dual5_exp(dual5 a)
{
    return dual5_apply(a, dual4_exp, dual4_exp_deriv1);
}

static inline dual5
dual5_log(dual5 a)
{
    return dual5_apply(a, dual4_log, dual4_log_deriv1);
}

/* Real analytic power: requires a positive real part for log(a). */
/* Integer powers by squaring; negative exponents invert the base once. */
static inline dual5
dual5_ipow(dual5 a, int64_t b)
{
    dual5 result = dual5_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return dual5_mul(a, a);
    }
    if (b == -1) {
        return dual5_inv(a);
    }
    if (b < 0) {
        a = dual5_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = dual5_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = dual5_mul(a, a);
        }
    }
    return result;
}

static inline dual5
dual5_rpow(dual5 a, dual5 b)
{
    return dual5_exp(dual5_mul(dual5_log(a), b));
}

static inline dual5
dual5_pow(dual5 a, dual5 b)
{
    return dual5_rpow(a, b);
}

static inline dual5
dual5_lpow(dual5 a, dual5 b)
{
    return dual5_rpow(a, b);
}

static inline dual5
dual5_sqrt(dual5 a)
{
    return dual5_apply(a, dual4_sqrt, dual4_sqrt_deriv1);
}

static inline dual5
dual5_cos(dual5 a)
{
    return dual5_apply(a, dual4_cos, dual4_cos_deriv1);
}

static inline dual5
dual5_sin(dual5 a)
{
    return dual5_apply(a, dual4_sin, dual4_sin_deriv1);
}

static inline dual5
dual5_tan(dual5 a)
{
    return dual5_apply(a, dual4_tan, dual4_tan_deriv1);
}

static inline dual5
dual5_acos(dual5 a)
{
    return dual5_apply(a, dual4_acos, dual4_acos_deriv1);
}

static inline dual5
dual5_asin(dual5 a)
{
    return dual5_apply(a, dual4_asin, dual4_asin_deriv1);
}

static inline dual5
dual5_atan(dual5 a)
{
    return dual5_apply(a, dual4_atan, dual4_atan_deriv1);
}

static inline dual5
dual5_cosh(dual5 a)
{
    return dual5_apply(a, dual4_cosh, dual4_cosh_deriv1);
}

static inline dual5
dual5_sinh(dual5 a)
{
    return dual5_apply(a, dual4_sinh, dual4_sinh_deriv1);
}

static inline dual5
dual5_tanh(dual5 a)
{
    return dual5_apply(a, dual4_tanh, dual4_tanh_deriv1);
}

static inline dual5
dual5_acosh(dual5 a)
{
    return dual5_apply(a, dual4_acosh, dual4_acosh_deriv1);
}

static inline dual5
dual5_asinh(dual5 a)
{
    return dual5_apply(a, dual4_asinh, dual4_asinh_deriv1);
}

static inline dual5
dual5_atanh(dual5 a)
{
    return dual5_apply(a, dual4_atanh, dual4_atanh_deriv1);
}

static inline dual5
dual5_expm1(dual5 a)
{
    return dual5_apply(a, dual4_expm1, dual4_expm1_deriv1);
}

static inline dual5
dual5_log1p(dual5 a)
{
    return dual5_apply(a, dual4_log1p, dual4_log1p_deriv1);
}

static inline dual5
dual5_radd(double a, dual5 b)
{
    return dual5_add(dual5_from_re(a), b);
}

static inline dual5
dual5_ladd(dual5 a, double b)
{
    return dual5_add(a, dual5_from_re(b));
}

static inline dual5
dual5_rsub(double a, dual5 b)
{
    return dual5_sub(dual5_from_re(a), b);
}

static inline dual5
dual5_lsub(dual5 a, double b)
{
    return dual5_sub(a, dual5_from_re(b));
}

static inline dual5
dual5_rmul(double a, dual5 b)
{
    return vec32_scale(b, a);
}

static inline dual5
dual5_lmul(dual5 a, double b)
{
    return vec32_scale(a, b);
}

static inline dual5
dual5_rdiv(double a, dual5 b)
{
    return dual5_div(dual5_from_re(a), b);
}

static inline dual5
dual5_ldiv(dual5 a, double b)
{
    return vec32_div_scalar(a, b);
}

static inline dual5
dual5_rpows(double a, dual5 b)
{
    return dual5_rpow(dual5_from_re(a), b);
}

static inline dual5
dual5_lpows(dual5 a, double b)
{
    return dual5_rpow(a, dual5_from_re(b));
}

static inline dual5
dual5_mix(dual5 a, dual5 b, double t)
{
    dual5 delta = dual5_mul(dual5_inv(a), b);
    return dual5_mul(a, dual5_lpows(delta, t));
}

static inline void
dual5_basis(dual5 result[ALGPY_DUAL5_N])
{
    int i, j;
    for (i = 0; i < ALGPY_DUAL5_N; i++) {
        for (j = 0; j < ALGPY_DUAL5_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
dual5_gens(dual5 result[5])
{
    result[0] = dual5_make((const double[ALGPY_DUAL5_N]){[1] = 1});
    result[1] = dual5_make((const double[ALGPY_DUAL5_N]){[2] = 1});
    result[2] = dual5_make((const double[ALGPY_DUAL5_N]){[4] = 1});
    result[3] = dual5_make((const double[ALGPY_DUAL5_N]){[8] = 1});
    result[4] = dual5_make((const double[ALGPY_DUAL5_N]){[16] = 1});
}

static inline void
dual5_coeffs(dual5 result[ALGPY_DUAL5_N*ALGPY_DUAL5_N])
{
    dual5 basis[ALGPY_DUAL5_N];
    int i, j;
    dual5_basis(basis);
    for (i = 0; i < ALGPY_DUAL5_N; i++) {
        for (j = 0; j < ALGPY_DUAL5_N; j++) {
            result[i*ALGPY_DUAL5_N+j] = dual5_mul(basis[i], basis[j]);
        }
    }
}

static inline void
dual5_minpoly(dual5 a, double result[33])
{
    double x = a.data[0];
    double x2 = x*x;
    double x4 = x2*x2;
    double x8 = x4*x4;
    double x16 = x8*x8;
    result[0] = x16*x16;
    result[1] = -32*x16*x8*x4*x2*x;
    result[2] = 496*x16*x8*x4*x2;
    result[3] = -4960*x16*x8*x4*x;
    result[4] = 35960*x16*x8*x4;
    result[5] = -201376*x16*x8*x2*x;
    result[6] = 906192*x16*x8*x2;
    result[7] = -3365856*x16*x8*x;
    result[8] = 10518300*x16*x8;
    result[9] = -28048800*x16*x4*x2*x;
    result[10] = 64512240*x16*x4*x2;
    result[11] = -129024480*x16*x4*x;
    result[12] = 225792840*x16*x4;
    result[13] = -347373600*x16*x2*x;
    result[14] = 471435600*x16*x2;
    result[15] = -565722720*x16*x;
    result[16] = 601080390*x16;
    result[17] = -565722720*x8*x4*x2*x;
    result[18] = 471435600*x8*x4*x2;
    result[19] = -347373600*x8*x4*x;
    result[20] = 225792840*x8*x4;
    result[21] = -129024480*x8*x2*x;
    result[22] = 64512240*x8*x2;
    result[23] = -28048800*x8*x;
    result[24] = 10518300*x8;
    result[25] = -3365856*x4*x2*x;
    result[26] = 906192*x4*x2;
    result[27] = -201376*x4*x;
    result[28] = 35960*x4;
    result[29] = -4960*x2*x;
    result[30] = 496*x2;
    result[31] = -32*x;
    result[32] = 1;
}

static inline void
dual5_eivals(dual5 a, complex result[32])
{
    int i;
    for (i = 0; i < ALGPY_DUAL5_N; i++) {
        result[i] = complex_from_re(a.data[0]);
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
dual5_lmat(dual5 a, double result[ALGPY_DUAL5_N*ALGPY_DUAL5_N])
{
    dual5 basis = {{0}};
    dual5 product;
    int row, col;
    for (col = 0; col < ALGPY_DUAL5_N; col++) {
        basis.data[col] = 1;
        product = dual5_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_DUAL5_N; row++) {
            result[row * ALGPY_DUAL5_N + col] = product.data[row];
        }
    }
}

static inline void
dual5_rmat(dual5 a, double result[ALGPY_DUAL5_N*ALGPY_DUAL5_N])
{
    dual5 basis = {{0}};
    dual5 product;
    int row, col;
    for (row = 0; row < ALGPY_DUAL5_N; row++) {
        basis.data[row] = 1;
        product = dual5_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_DUAL5_N; col++) {
            result[row * ALGPY_DUAL5_N + col] = product.data[col];
        }
    }
}



/* Additional mixed derivatives, kept for future extensions. */
static inline dual5
dual5_exp_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_exp_deriv1, dual4_exp_deriv2);
}

static inline dual5
dual5_expm1_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_expm1_deriv1, dual4_expm1_deriv2);
}

static inline dual5
dual5_log_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_log_deriv1, dual4_log_deriv2);
}

static inline dual5
dual5_log1p_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_log1p_deriv1, dual4_log1p_deriv2);
}

static inline dual5
dual5_inv_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_inv_deriv1, dual4_inv_deriv2);
}

static inline dual5
dual5_square_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_square_deriv1, dual4_square_deriv2);
}

static inline dual5
dual5_sqrt_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_sqrt_deriv1, dual4_sqrt_deriv2);
}

static inline dual5
dual5_sin_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_sin_deriv1, dual4_sin_deriv2);
}

static inline dual5
dual5_cos_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_cos_deriv1, dual4_cos_deriv2);
}

static inline dual5
dual5_tan_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_tan_deriv1, dual4_tan_deriv2);
}

static inline dual5
dual5_asin_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_asin_deriv1, dual4_asin_deriv2);
}

static inline dual5
dual5_acos_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_acos_deriv1, dual4_acos_deriv2);
}

static inline dual5
dual5_atan_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_atan_deriv1, dual4_atan_deriv2);
}

static inline dual5
dual5_sinh_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_sinh_deriv1, dual4_sinh_deriv2);
}

static inline dual5
dual5_cosh_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_cosh_deriv1, dual4_cosh_deriv2);
}

static inline dual5
dual5_tanh_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_tanh_deriv1, dual4_tanh_deriv2);
}

static inline dual5
dual5_asinh_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_asinh_deriv1, dual4_asinh_deriv2);
}

static inline dual5
dual5_acosh_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_acosh_deriv1, dual4_acosh_deriv2);
}

static inline dual5
dual5_atanh_deriv1(dual5 a)
{
    return dual5_apply(a, dual4_atanh_deriv1, dual4_atanh_deriv2);
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline dual5
dual5_grade0(dual5 a)
{
    const double *x = a.data;
    return dual5_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline dual5
dual5_grade1(dual5 a)
{
    const double *x = a.data;
    return dual5_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline dual5
dual5_grade2(dual5 a)
{
    const double *x = a.data;
    return dual5_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline dual5
dual5_grade3(dual5 a)
{
    const double *x = a.data;
    return dual5_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline dual5
dual5_grade4(dual5 a)
{
    const double *x = a.data;
    return dual5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline dual5
dual5_grade5(dual5 a)
{
    const double *x = a.data;
    return dual5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

#endif /* ALGPY_DUAL5_H */
