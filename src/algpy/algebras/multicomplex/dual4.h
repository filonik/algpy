#ifndef ALGPY_DUAL4_H
#define ALGPY_DUAL4_H

#include <string.h>
#include <stdint.h>
#include "../../vectorspaces/vec16.h"
#include "dual3.h"

#define ALGPY_DUAL4_N ALGPY_VEC16_N

typedef vec16 dual4;

static inline dual4
dual4_make(const double coeffs[ALGPY_DUAL4_N])
{
    dual4 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

/* For a = u + v*e with e*e = 0, f(a) = f(u) + v*f'(u)*e.
 * Both callbacks are lower-dimensional dual functions, carrying all mixed
 * derivatives. They must be defined and differentiable at the real part. */
static inline dual4
dual4_apply(dual4 a, dual3 (*f)(dual3), dual3 (*deriv1)(dual3))
{
    dual3 u = dual3_make(a.data);
    dual3 v = dual3_make(a.data + 8);
    dual3 value = f(u);
    dual3 tangent = dual3_mul(v, deriv1(u));
    dual4 result;
    memcpy(result.data, value.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_DUAL4_N/2, tangent.data, sizeof(result.data)/2);
    return result;
}

static inline dual4
dual4_from_re(double a)
{
    return dual4_make((const double[ALGPY_DUAL4_N]){a});
}

static inline double
dual4_to_re(dual4 a)
{
    return a.data[0];
}

static inline dual4
dual4_zero(void)
{
    return vec16_zero();
}

static inline dual4
dual4_one(void)
{
    return dual4_from_re(1);
}

static inline int
dual4_is_real(dual4 a)
{
    return vec16_is_real(a);
}

static inline dual4
dual4_pos(dual4 a)
{
    return vec16_pos(a);
}

static inline dual4
dual4_neg(dual4 a)
{
    return vec16_neg(a);
}

static inline dual4
dual4_re(dual4 a)
{
    return dual4_from_re(a.data[0]);
}

static inline dual4
dual4_im(dual4 a)
{
    return vec16_sub(a, dual4_re(a));
}

static inline dual4
dual4_conjugate(dual4 a)
{
    return dual4_make((const double[]){a.data[0], a.data[1], a.data[2], a.data[3], a.data[4], a.data[5], a.data[6], a.data[7], -a.data[8], -a.data[9], -a.data[10], -a.data[11], -a.data[12], -a.data[13], -a.data[14], -a.data[15]});
}

static inline dual4
dual4_add(dual4 a, dual4 b)
{
    return vec16_add(a, b);
}

static inline dual4
dual4_sub(dual4 a, dual4 b)
{
    return vec16_sub(a, b);
}

static inline int
dual4_equal(dual4 a, dual4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
dual4_not_equal(dual4 a, dual4 b)
{
    return !dual4_equal(a, b);
}

static inline dual4
dual4_mul(dual4 a, dual4 b)
{
    if (dual4_is_real(b)) {
        return vec16_scale(a, b.data[0]);
    }
    if (dual4_is_real(a)) {
        return vec16_scale(b, a.data[0]);
    }
    dual3 u = dual3_make(a.data);
    dual3 v = dual3_make(a.data + 8);
    dual3 s = dual3_make(b.data);
    dual3 t = dual3_make(b.data + 8);
    dual3 lo = dual3_mul(u, s);
    dual3 hi = dual3_add(dual3_mul(u, t), dual3_mul(v, s));
    dual4 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_DUAL4_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

/* Norm a*conjugate(a) = u*u, since the outer generator squares to zero. */
static inline dual3
dual4_norm2(dual4 a)
{
    dual3 u = dual3_make(a.data);
    return dual3_square(u);
}

static inline dual2
dual4_norm4(dual4 a)
{
    return dual3_norm2(dual4_norm2(a));
}

static inline dual1
dual4_norm8(dual4 a)
{
    return dual2_norm2(dual4_norm4(a));
}

static inline double
dual4_norm16(dual4 a)
{
    return dual1_norm2(dual4_norm8(a));
}

static inline double
dual4_norm(dual4 a)
{
    return fabs(a.data[0]);
}

static inline dual4
dual4_inv(dual4 a)
{
    dual3 u = dual3_make(a.data);
    dual3 v = dual3_make(a.data + 8);
    dual3 lo = dual3_inv(u);
    dual3 hi = dual3_neg(dual3_mul(dual3_mul(v, lo), lo));
    dual4 result;
    memcpy(result.data, lo.data, sizeof(result.data)/2);
    memcpy(result.data + ALGPY_DUAL4_N/2, hi.data, sizeof(result.data)/2);
    return result;
}

static inline dual4
dual4_div(dual4 a, dual4 b)
{
    if (dual4_is_real(b)) {
        return vec16_div_scalar(a, b.data[0]);
    }
    return dual4_mul(a, dual4_inv(b));
}

static inline dual4
dual4_square(dual4 a)
{
    return dual4_mul(a, a);
}

static inline dual4
dual4_exp(dual4 a)
{
    return dual4_apply(a, dual3_exp, dual3_exp_deriv1);
}

static inline dual4
dual4_log(dual4 a)
{
    return dual4_apply(a, dual3_log, dual3_log_deriv1);
}

/* Real analytic power: requires a positive real part for log(a). */
/* Integer powers by squaring; negative exponents invert the base once. */
static inline dual4
dual4_ipow(dual4 a, int64_t b)
{
    dual4 result = dual4_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return dual4_mul(a, a);
    }
    if (b == -1) {
        return dual4_inv(a);
    }
    if (b < 0) {
        a = dual4_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = dual4_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = dual4_mul(a, a);
        }
    }
    return result;
}

static inline dual4
dual4_rpow(dual4 a, dual4 b)
{
    return dual4_exp(dual4_mul(dual4_log(a), b));
}

static inline dual4
dual4_pow(dual4 a, dual4 b)
{
    return dual4_rpow(a, b);
}

static inline dual4
dual4_lpow(dual4 a, dual4 b)
{
    return dual4_rpow(a, b);
}

static inline dual4
dual4_sqrt(dual4 a)
{
    return dual4_apply(a, dual3_sqrt, dual3_sqrt_deriv1);
}

static inline dual4
dual4_cos(dual4 a)
{
    return dual4_apply(a, dual3_cos, dual3_cos_deriv1);
}

static inline dual4
dual4_sin(dual4 a)
{
    return dual4_apply(a, dual3_sin, dual3_sin_deriv1);
}

static inline dual4
dual4_tan(dual4 a)
{
    return dual4_apply(a, dual3_tan, dual3_tan_deriv1);
}

static inline dual4
dual4_acos(dual4 a)
{
    return dual4_apply(a, dual3_acos, dual3_acos_deriv1);
}

static inline dual4
dual4_asin(dual4 a)
{
    return dual4_apply(a, dual3_asin, dual3_asin_deriv1);
}

static inline dual4
dual4_atan(dual4 a)
{
    return dual4_apply(a, dual3_atan, dual3_atan_deriv1);
}

static inline dual4
dual4_cosh(dual4 a)
{
    return dual4_apply(a, dual3_cosh, dual3_cosh_deriv1);
}

static inline dual4
dual4_sinh(dual4 a)
{
    return dual4_apply(a, dual3_sinh, dual3_sinh_deriv1);
}

static inline dual4
dual4_tanh(dual4 a)
{
    return dual4_apply(a, dual3_tanh, dual3_tanh_deriv1);
}

static inline dual4
dual4_acosh(dual4 a)
{
    return dual4_apply(a, dual3_acosh, dual3_acosh_deriv1);
}

static inline dual4
dual4_asinh(dual4 a)
{
    return dual4_apply(a, dual3_asinh, dual3_asinh_deriv1);
}

static inline dual4
dual4_atanh(dual4 a)
{
    return dual4_apply(a, dual3_atanh, dual3_atanh_deriv1);
}

static inline dual4
dual4_expm1(dual4 a)
{
    return dual4_apply(a, dual3_expm1, dual3_expm1_deriv1);
}

static inline dual4
dual4_log1p(dual4 a)
{
    return dual4_apply(a, dual3_log1p, dual3_log1p_deriv1);
}

static inline dual4
dual4_radd(double a, dual4 b)
{
    return dual4_add(dual4_from_re(a), b);
}

static inline dual4
dual4_ladd(dual4 a, double b)
{
    return dual4_add(a, dual4_from_re(b));
}

static inline dual4
dual4_rsub(double a, dual4 b)
{
    return dual4_sub(dual4_from_re(a), b);
}

static inline dual4
dual4_lsub(dual4 a, double b)
{
    return dual4_sub(a, dual4_from_re(b));
}

static inline dual4
dual4_rmul(double a, dual4 b)
{
    return vec16_scale(b, a);
}

static inline dual4
dual4_lmul(dual4 a, double b)
{
    return vec16_scale(a, b);
}

static inline dual4
dual4_rdiv(double a, dual4 b)
{
    return dual4_div(dual4_from_re(a), b);
}

static inline dual4
dual4_ldiv(dual4 a, double b)
{
    return vec16_div_scalar(a, b);
}

static inline dual4
dual4_rpows(double a, dual4 b)
{
    return dual4_rpow(dual4_from_re(a), b);
}

static inline dual4
dual4_lpows(dual4 a, double b)
{
    return dual4_rpow(a, dual4_from_re(b));
}

static inline dual4
dual4_mix(dual4 a, dual4 b, double t)
{
    dual4 delta = dual4_mul(dual4_inv(a), b);
    return dual4_mul(a, dual4_lpows(delta, t));
}

static inline void
dual4_basis(dual4 result[ALGPY_DUAL4_N])
{
    int i, j;
    for (i = 0; i < ALGPY_DUAL4_N; i++) {
        for (j = 0; j < ALGPY_DUAL4_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
dual4_gens(dual4 result[4])
{
    result[0] = dual4_make((const double[ALGPY_DUAL4_N]){[1] = 1});
    result[1] = dual4_make((const double[ALGPY_DUAL4_N]){[2] = 1});
    result[2] = dual4_make((const double[ALGPY_DUAL4_N]){[4] = 1});
    result[3] = dual4_make((const double[ALGPY_DUAL4_N]){[8] = 1});
}

static inline void
dual4_coeffs(dual4 result[ALGPY_DUAL4_N*ALGPY_DUAL4_N])
{
    dual4 basis[ALGPY_DUAL4_N];
    int i, j;
    dual4_basis(basis);
    for (i = 0; i < ALGPY_DUAL4_N; i++) {
        for (j = 0; j < ALGPY_DUAL4_N; j++) {
            result[i*ALGPY_DUAL4_N+j] = dual4_mul(basis[i], basis[j]);
        }
    }
}

static inline void
dual4_minpoly(dual4 a, double result[17])
{
    double x = a.data[0];
    double x2 = x*x;
    double x4 = x2*x2;
    double x8 = x4*x4;
    result[0] = x8*x8;
    result[1] = -16*x8*x4*x2*x;
    result[2] = 120*x8*x4*x2;
    result[3] = -560*x8*x4*x;
    result[4] = 1820*x8*x4;
    result[5] = -4368*x8*x2*x;
    result[6] = 8008*x8*x2;
    result[7] = -11440*x8*x;
    result[8] = 12870*x8;
    result[9] = -11440*x4*x2*x;
    result[10] = 8008*x4*x2;
    result[11] = -4368*x4*x;
    result[12] = 1820*x4;
    result[13] = -560*x2*x;
    result[14] = 120*x2;
    result[15] = -16*x;
    result[16] = 1;
}

static inline void
dual4_eivals(dual4 a, complex result[16])
{
    int i;
    for (i = 0; i < ALGPY_DUAL4_N; i++) {
        result[i] = complex_from_re(a.data[0]);
    }
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
dual4_lmat(dual4 a, double result[ALGPY_DUAL4_N*ALGPY_DUAL4_N])
{
    dual4 basis = {{0}};
    dual4 product;
    int row, col;
    for (col = 0; col < ALGPY_DUAL4_N; col++) {
        basis.data[col] = 1;
        product = dual4_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_DUAL4_N; row++) {
            result[row * ALGPY_DUAL4_N + col] = product.data[row];
        }
    }
}

static inline void
dual4_rmat(dual4 a, double result[ALGPY_DUAL4_N*ALGPY_DUAL4_N])
{
    dual4 basis = {{0}};
    dual4 product;
    int row, col;
    for (row = 0; row < ALGPY_DUAL4_N; row++) {
        basis.data[row] = 1;
        product = dual4_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_DUAL4_N; col++) {
            result[row * ALGPY_DUAL4_N + col] = product.data[col];
        }
    }
}


/* Additional mixed derivatives for dual5. */
static inline dual4
dual4_exp_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_exp_deriv1, dual3_exp_deriv2);
}

static inline dual4
dual4_expm1_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_expm1_deriv1, dual3_expm1_deriv2);
}

static inline dual4
dual4_log_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_log_deriv1, dual3_log_deriv2);
}

static inline dual4
dual4_log1p_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_log1p_deriv1, dual3_log1p_deriv2);
}

static inline dual4
dual4_inv_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_inv_deriv1, dual3_inv_deriv2);
}

static inline dual4
dual4_square_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_square_deriv1, dual3_square_deriv2);
}

static inline dual4
dual4_sqrt_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_sqrt_deriv1, dual3_sqrt_deriv2);
}

static inline dual4
dual4_sin_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_sin_deriv1, dual3_sin_deriv2);
}

static inline dual4
dual4_cos_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_cos_deriv1, dual3_cos_deriv2);
}

static inline dual4
dual4_tan_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_tan_deriv1, dual3_tan_deriv2);
}

static inline dual4
dual4_asin_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_asin_deriv1, dual3_asin_deriv2);
}

static inline dual4
dual4_acos_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_acos_deriv1, dual3_acos_deriv2);
}

static inline dual4
dual4_atan_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_atan_deriv1, dual3_atan_deriv2);
}

static inline dual4
dual4_sinh_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_sinh_deriv1, dual3_sinh_deriv2);
}

static inline dual4
dual4_cosh_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_cosh_deriv1, dual3_cosh_deriv2);
}

static inline dual4
dual4_tanh_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_tanh_deriv1, dual3_tanh_deriv2);
}

static inline dual4
dual4_asinh_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_asinh_deriv1, dual3_asinh_deriv2);
}

static inline dual4
dual4_acosh_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_acosh_deriv1, dual3_acosh_deriv2);
}

static inline dual4
dual4_atanh_deriv1(dual4 a)
{
    return dual4_apply(a, dual3_atanh_deriv1, dual3_atanh_deriv2);
}


/* Additional mixed derivatives, kept for future extensions. */
static inline dual4
dual4_exp_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_exp_deriv2, dual3_exp_deriv3);
}

static inline dual4
dual4_expm1_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_expm1_deriv2, dual3_expm1_deriv3);
}

static inline dual4
dual4_log_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_log_deriv2, dual3_log_deriv3);
}

static inline dual4
dual4_log1p_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_log1p_deriv2, dual3_log1p_deriv3);
}

static inline dual4
dual4_inv_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_inv_deriv2, dual3_inv_deriv3);
}

static inline dual4
dual4_square_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_square_deriv2, dual3_square_deriv3);
}

static inline dual4
dual4_sqrt_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_sqrt_deriv2, dual3_sqrt_deriv3);
}

static inline dual4
dual4_sin_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_sin_deriv2, dual3_sin_deriv3);
}

static inline dual4
dual4_cos_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_cos_deriv2, dual3_cos_deriv3);
}

static inline dual4
dual4_tan_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_tan_deriv2, dual3_tan_deriv3);
}

static inline dual4
dual4_asin_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_asin_deriv2, dual3_asin_deriv3);
}

static inline dual4
dual4_acos_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_acos_deriv2, dual3_acos_deriv3);
}

static inline dual4
dual4_atan_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_atan_deriv2, dual3_atan_deriv3);
}

static inline dual4
dual4_sinh_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_sinh_deriv2, dual3_sinh_deriv3);
}

static inline dual4
dual4_cosh_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_cosh_deriv2, dual3_cosh_deriv3);
}

static inline dual4
dual4_tanh_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_tanh_deriv2, dual3_tanh_deriv3);
}

static inline dual4
dual4_asinh_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_asinh_deriv2, dual3_asinh_deriv3);
}

static inline dual4
dual4_acosh_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_acosh_deriv2, dual3_acosh_deriv3);
}

static inline dual4
dual4_atanh_deriv2(dual4 a)
{
    return dual4_apply(a, dual3_atanh_deriv2, dual3_atanh_deriv3);
}

/* Grade projections: grade k keeps the coefficients whose basis index has
 * k generator bits. */
static inline dual4
dual4_grade0(dual4 a)
{
    const double *x = a.data;
    return dual4_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline dual4
dual4_grade1(dual4 a)
{
    const double *x = a.data;
    return dual4_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0});
}

static inline dual4
dual4_grade2(dual4 a)
{
    const double *x = a.data;
    return dual4_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0});
}

static inline dual4
dual4_grade3(dual4 a)
{
    const double *x = a.data;
    return dual4_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0});
}

static inline dual4
dual4_grade4(dual4 a)
{
    const double *x = a.data;
    return dual4_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15]});
}

#endif /* ALGPY_DUAL4_H */
