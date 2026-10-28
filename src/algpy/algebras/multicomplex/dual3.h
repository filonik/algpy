#ifndef ALGPY_DUAL3_H
#define ALGPY_DUAL3_H

#include "../../vectorspaces/vec8.h"
#include "dual2.h"

#define ALGPY_DUAL3_N ALGPY_VEC8_N

typedef vec8 dual3;

static inline dual3
dual3_make(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
    dual3 result = {{x0, x1, x2, x3, x4, x5, x6, x7}};
    return result;
}

/* For a = u + v*e with e*e = 0, f(a) = f(u) + v*f'(u)*e.
 * Both callbacks are lower-dimensional dual functions, carrying all mixed
 * derivatives. They must be defined and differentiable at the real part. */
static inline dual3
dual3_apply(dual3 a, dual2 (*f)(dual2), dual2 (*deriv1)(dual2))
{
    dual2 u = dual2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    dual2 v = dual2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    dual2 value = f(u);
    dual2 tangent = dual2_mul(v, deriv1(u));
    return dual3_make(value.data[0], value.data[1], value.data[2], value.data[3], tangent.data[0], tangent.data[1], tangent.data[2], tangent.data[3]);
}

static inline dual3
dual3_from_re(double a)
{
    return dual3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
dual3_to_re(dual3 a)
{
    return a.data[0];
}

static inline dual3
dual3_zero(void)
{
    return vec8_zero();
}

static inline dual3
dual3_one(void)
{
    return dual3_from_re(1);
}

static inline int
dual3_is_real(dual3 a)
{
    return vec8_is_real(a);
}

static inline dual3
dual3_pos(dual3 a)
{
    return vec8_pos(a);
}

static inline dual3
dual3_neg(dual3 a)
{
    return vec8_neg(a);
}

static inline dual3
dual3_re(dual3 a)
{
    return dual3_from_re(a.data[0]);
}

static inline dual3
dual3_im(dual3 a)
{
    return vec8_sub(a, dual3_re(a));
}

static inline dual3
dual3_conjugate(dual3 a)
{
    return dual3_make(a.data[0], a.data[1], a.data[2], a.data[3], -a.data[4], -a.data[5], -a.data[6], -a.data[7]);
}

static inline dual3
dual3_add(dual3 a, dual3 b)
{
    return vec8_add(a, b);
}

static inline dual3
dual3_sub(dual3 a, dual3 b)
{
    return vec8_sub(a, b);
}

static inline int
dual3_equal(dual3 a, dual3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
dual3_not_equal(dual3 a, dual3 b)
{
    return !dual3_equal(a, b);
}

static inline dual3
dual3_mul(dual3 a, dual3 b)
{
    if (dual3_is_real(b)) {
        return vec8_scale(a, b.data[0]);
    }
    if (dual3_is_real(a)) {
        return vec8_scale(b, a.data[0]);
    }
    dual2 u = dual2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    dual2 v = dual2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    dual2 s = dual2_make(b.data[0], b.data[1], b.data[2], b.data[3]);
    dual2 t = dual2_make(b.data[4], b.data[5], b.data[6], b.data[7]);
    dual2 lo = dual2_mul(u, s);
    dual2 hi = dual2_add(dual2_mul(u, t), dual2_mul(v, s));
    return dual3_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], hi.data[0], hi.data[1], hi.data[2], hi.data[3]);
}

/* Norm a*conjugate(a) = u*u, since the outer generator squares to zero. */
static inline dual2
dual3_norm2(dual3 a)
{
    dual2 u = dual2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    return dual2_square(u);
}

static inline dual1
dual3_norm4(dual3 a)
{
    return dual2_norm2(dual3_norm2(a));
}

static inline double
dual3_norm8(dual3 a)
{
    return dual1_norm2(dual3_norm4(a));
}

static inline double
dual3_norm(dual3 a)
{
    return fabs(a.data[0]);
}

static inline dual3
dual3_inv(dual3 a)
{
    dual2 u = dual2_make(a.data[0], a.data[1], a.data[2], a.data[3]);
    dual2 v = dual2_make(a.data[4], a.data[5], a.data[6], a.data[7]);
    dual2 lo = dual2_inv(u);
    dual2 hi = dual2_neg(dual2_mul(dual2_mul(v, lo), lo));
    return dual3_make(lo.data[0], lo.data[1], lo.data[2], lo.data[3], hi.data[0], hi.data[1], hi.data[2], hi.data[3]);
}

static inline dual3
dual3_div(dual3 a, dual3 b)
{
    if (dual3_is_real(b)) {
        return vec8_div_scalar(a, b.data[0]);
    }
    return dual3_mul(a, dual3_inv(b));
}

static inline dual3
dual3_square(dual3 a)
{
    return dual3_mul(a, a);
}

static inline dual3
dual3_exp(dual3 a)
{
    return dual3_apply(a, dual2_exp, dual2_exp_deriv1);
}

static inline dual3
dual3_log(dual3 a)
{
    return dual3_apply(a, dual2_log, dual2_log_deriv1);
}

/* Real analytic power: requires a positive real part for log(a). */
static inline dual3
dual3_pow(dual3 a, dual3 b)
{
    return dual3_exp(dual3_mul(dual3_log(a), b));
}

static inline dual3
dual3_sqrt(dual3 a)
{
    return dual3_apply(a, dual2_sqrt, dual2_sqrt_deriv1);
}

static inline dual3
dual3_cos(dual3 a)
{
    return dual3_apply(a, dual2_cos, dual2_cos_deriv1);
}

static inline dual3
dual3_sin(dual3 a)
{
    return dual3_apply(a, dual2_sin, dual2_sin_deriv1);
}

static inline dual3
dual3_tan(dual3 a)
{
    return dual3_apply(a, dual2_tan, dual2_tan_deriv1);
}

static inline dual3
dual3_acos(dual3 a)
{
    return dual3_apply(a, dual2_acos, dual2_acos_deriv1);
}

static inline dual3
dual3_asin(dual3 a)
{
    return dual3_apply(a, dual2_asin, dual2_asin_deriv1);
}

static inline dual3
dual3_atan(dual3 a)
{
    return dual3_apply(a, dual2_atan, dual2_atan_deriv1);
}

static inline dual3
dual3_cosh(dual3 a)
{
    return dual3_apply(a, dual2_cosh, dual2_cosh_deriv1);
}

static inline dual3
dual3_sinh(dual3 a)
{
    return dual3_apply(a, dual2_sinh, dual2_sinh_deriv1);
}

static inline dual3
dual3_tanh(dual3 a)
{
    return dual3_apply(a, dual2_tanh, dual2_tanh_deriv1);
}

static inline dual3
dual3_acosh(dual3 a)
{
    return dual3_apply(a, dual2_acosh, dual2_acosh_deriv1);
}

static inline dual3
dual3_asinh(dual3 a)
{
    return dual3_apply(a, dual2_asinh, dual2_asinh_deriv1);
}

static inline dual3
dual3_atanh(dual3 a)
{
    return dual3_apply(a, dual2_atanh, dual2_atanh_deriv1);
}

static inline dual3
dual3_expm1(dual3 a)
{
    return dual3_apply(a, dual2_expm1, dual2_expm1_deriv1);
}

static inline dual3
dual3_log1p(dual3 a)
{
    return dual3_apply(a, dual2_log1p, dual2_log1p_deriv1);
}

static inline dual3
dual3_radd(double a, dual3 b)
{
    return dual3_add(dual3_from_re(a), b);
}

static inline dual3
dual3_ladd(dual3 a, double b)
{
    return dual3_add(a, dual3_from_re(b));
}

static inline dual3
dual3_rsub(double a, dual3 b)
{
    return dual3_sub(dual3_from_re(a), b);
}

static inline dual3
dual3_lsub(dual3 a, double b)
{
    return dual3_sub(a, dual3_from_re(b));
}

static inline dual3
dual3_rmul(double a, dual3 b)
{
    return vec8_scale(b, a);
}

static inline dual3
dual3_lmul(dual3 a, double b)
{
    return vec8_scale(a, b);
}

static inline dual3
dual3_rdiv(double a, dual3 b)
{
    return dual3_div(dual3_from_re(a), b);
}

static inline dual3
dual3_ldiv(dual3 a, double b)
{
    return vec8_div_scalar(a, b);
}

static inline dual3
dual3_rpow(double a, dual3 b)
{
    return dual3_pow(dual3_from_re(a), b);
}

static inline dual3
dual3_lpow(dual3 a, double b)
{
    return dual3_pow(a, dual3_from_re(b));
}

static inline dual3
dual3_mix(dual3 a, dual3 b, double t)
{
    dual3 delta = dual3_mul(dual3_inv(a), b);
    return dual3_mul(a, dual3_lpow(delta, t));
}

static inline void
dual3_basis(dual3 result[ALGPY_DUAL3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_DUAL3_N; i++) {
        for (j = 0; j < ALGPY_DUAL3_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
dual3_coeffs(dual3 result[ALGPY_DUAL3_N*ALGPY_DUAL3_N])
{
    dual3 basis[ALGPY_DUAL3_N];
    int i, j;
    dual3_basis(basis);
    for (i = 0; i < ALGPY_DUAL3_N; i++) {
        for (j = 0; j < ALGPY_DUAL3_N; j++) {
            result[i*ALGPY_DUAL3_N+j] = dual3_mul(basis[i], basis[j]);
        }
    }
}

static inline void
dual3_minpoly(dual3 a, double result[9])
{
    double x = a.data[0];
    double x2 = x*x;
    double x4 = x2*x2;
    result[0] = x4*x4;
    result[1] = -8*x4*x2*x;
    result[2] = 28*x4*x2;
    result[3] = -56*x4*x;
    result[4] = 70*x4;
    result[5] = -56*x2*x;
    result[6] = 28*x2;
    result[7] = -8*x;
    result[8] = 1;
}

static inline void
dual3_eivals(dual3 a, complex result[8])
{
    int i;
    for (i = 0; i < ALGPY_DUAL3_N; i++) {
        result[i] = complex_from_re(a.data[0]);
    }
}

/* Derivatives through order 1; each dual level uses one more field derivative. */
static inline dual3
dual3_exp_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_exp_deriv1, dual2_exp_deriv2);
}

static inline dual3
dual3_expm1_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_expm1_deriv1, dual2_expm1_deriv2);
}

static inline dual3
dual3_log_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_log_deriv1, dual2_log_deriv2);
}

static inline dual3
dual3_log1p_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_log1p_deriv1, dual2_log1p_deriv2);
}

static inline dual3
dual3_inv_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_inv_deriv1, dual2_inv_deriv2);
}

static inline dual3
dual3_square_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_square_deriv1, dual2_square_deriv2);
}

static inline dual3
dual3_sqrt_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_sqrt_deriv1, dual2_sqrt_deriv2);
}

static inline dual3
dual3_sin_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_sin_deriv1, dual2_sin_deriv2);
}

static inline dual3
dual3_cos_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_cos_deriv1, dual2_cos_deriv2);
}

static inline dual3
dual3_tan_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_tan_deriv1, dual2_tan_deriv2);
}

static inline dual3
dual3_asin_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_asin_deriv1, dual2_asin_deriv2);
}

static inline dual3
dual3_acos_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_acos_deriv1, dual2_acos_deriv2);
}

static inline dual3
dual3_atan_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_atan_deriv1, dual2_atan_deriv2);
}

static inline dual3
dual3_sinh_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_sinh_deriv1, dual2_sinh_deriv2);
}

static inline dual3
dual3_cosh_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_cosh_deriv1, dual2_cosh_deriv2);
}

static inline dual3
dual3_tanh_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_tanh_deriv1, dual2_tanh_deriv2);
}

static inline dual3
dual3_asinh_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_asinh_deriv1, dual2_asinh_deriv2);
}

static inline dual3
dual3_acosh_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_acosh_deriv1, dual2_acosh_deriv2);
}

static inline dual3
dual3_atanh_deriv1(dual3 a)
{
    return dual3_apply(a, dual2_atanh_deriv1, dual2_atanh_deriv2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
dual3_lmat(dual3 a, double result[ALGPY_DUAL3_N*ALGPY_DUAL3_N])
{
    dual3 coeffs[ALGPY_DUAL3_N * ALGPY_DUAL3_N];
    int row, col, i;
    dual3_coeffs(coeffs);
    for (row = 0; row < ALGPY_DUAL3_N; row++) {
        for (col = 0; col < ALGPY_DUAL3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_DUAL3_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_DUAL3_N + col].data[row];
            }
            result[row * ALGPY_DUAL3_N + col] = value;
        }
    }
}

static inline void
dual3_rmat(dual3 a, double result[ALGPY_DUAL3_N*ALGPY_DUAL3_N])
{
    dual3 coeffs[ALGPY_DUAL3_N * ALGPY_DUAL3_N];
    int row, col, i;
    dual3_coeffs(coeffs);
    for (row = 0; row < ALGPY_DUAL3_N; row++) {
        for (col = 0; col < ALGPY_DUAL3_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_DUAL3_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_DUAL3_N + i].data[col];
            }
            result[row * ALGPY_DUAL3_N + col] = value;
        }
    }
}


#endif /* ALGPY_DUAL3_H */
