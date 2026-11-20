/* Projective geometric algebra in 2D.
 * Generator squares: 0, 1, 1.
 * 3 generators; 8 real coefficients. */
#ifndef ALGPY_PGA2_H
#define ALGPY_PGA2_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../polynomials/hermite.h"
#include "../../vectorspaces/vec8.h"
#include "../../lattices/clifford/exterior3.h"
#include "../multicomplex/dual1.h"

#define ALGPY_PGA2_N 8

typedef vec8 pga2;

static inline pga2
pga2_make(const double coeffs[ALGPY_PGA2_N])
{
    pga2 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline pga2
pga2_grade0(pga2 a)
{
    const double *x = a.data;
    return pga2_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0});
}

static inline pga2
pga2_grade1(pga2 a)
{
    const double *x = a.data;
    return pga2_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0});
}

static inline pga2
pga2_grade2(pga2 a)
{
    const double *x = a.data;
    return pga2_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0});
}

static inline pga2
pga2_grade3(pga2 a)
{
    const double *x = a.data;
    return pga2_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7]});
}

static inline pga2
pga2_from_re(double a)
{
    return pga2_make((const double[ALGPY_PGA2_N]){a});
}

static inline double
pga2_to_re(pga2 a)
{
    return a.data[0];
}

static inline pga2
pga2_zero(void)
{
    return vec8_zero();
}

static inline pga2
pga2_one(void)
{
    return pga2_from_re(1);
}

static inline pga2
pga2_pos(pga2 a)
{
    return vec8_pos(a);
}

static inline pga2
pga2_neg(pga2 a)
{
    return vec8_neg(a);
}

static inline pga2
pga2_add(pga2 a, pga2 b)
{
    return vec8_add(a, b);
}

static inline pga2
pga2_sub(pga2 a, pga2 b)
{
    return vec8_sub(a, b);
}

static inline int
pga2_equal(pga2 a, pga2 b)
{
    return vec8_all_eq(a, b);
}

static inline int
pga2_not_equal(pga2 a, pga2 b)
{
    return !pga2_equal(a, b);
}

static inline int
pga2_is_real(pga2 a)
{
    return vec8_is_real(a);
}

static inline pga2
pga2_re(pga2 a)
{
    return pga2_from_re(a.data[0]);
}

static inline pga2
pga2_im(pga2 a)
{
    return pga2_sub(a, pga2_re(a));
}

static inline pga2
pga2_mul(pga2 a, pga2 b)
{
    const double *x = a.data, *y = b.data;
    return pga2_make((const double[]){
        x[0]*y[0] + x[2]*y[2] + x[4]*y[4] - x[6]*y[6],
        x[0]*y[1] + x[1]*y[0] - x[2]*y[3] + x[3]*y[2] - x[4]*y[5] + x[5]*y[4] - x[6]*y[7] - x[7]*y[6],
        x[0]*y[2] + x[2]*y[0] - x[4]*y[6] + x[6]*y[4],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0] + x[4]*y[7] - x[5]*y[6] + x[6]*y[5] + x[7]*y[4],
        x[0]*y[4] + x[2]*y[6] + x[4]*y[0] - x[6]*y[2],
        x[0]*y[5] + x[1]*y[4] - x[2]*y[7] + x[3]*y[6] - x[4]*y[1] + x[5]*y[0] - x[6]*y[3] - x[7]*y[2],
        x[0]*y[6] + x[2]*y[4] - x[4]*y[2] + x[6]*y[0],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0]});
}

static inline pga2
pga2_square(pga2 a)
{
    const double *x = a.data;
    return pga2_make((const double[]){
        x[0]*x[0] + x[2]*x[2] + x[4]*x[4] - x[6]*x[6],
        2*x[0]*x[1] - 2*x[6]*x[7],
        2*x[0]*x[2],
        2*x[0]*x[3] + 2*x[4]*x[7],
        2*x[0]*x[4],
        2*x[0]*x[5] - 2*x[2]*x[7],
        2*x[0]*x[6],
        2*x[0]*x[7] + 2*x[1]*x[6] - 2*x[2]*x[5] + 2*x[3]*x[4]});
}

static inline pga2
pga2_conjugate(pga2 a)
{
    const double *x = a.data;
    return pga2_make((const double[]){x[0], -x[1], -x[2], -x[3], -x[4], -x[5], -x[6], x[7]});
}

static inline void
pga2_basis(pga2 result[8])
{
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
pga2_gens(pga2 result[3])
{
    result[0] = pga2_make((const double[ALGPY_PGA2_N]){[1] = 1});
    result[1] = pga2_make((const double[ALGPY_PGA2_N]){[2] = 1});
    result[2] = pga2_make((const double[ALGPY_PGA2_N]){[4] = 1});
}

static inline void
pga2_coeffs(pga2 result[64])
{
    pga2 basis[8];
    int i, j;
    pga2_basis(basis);
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            result[i*8+j] = pga2_mul(basis[i], basis[j]);
        }
    }
}

static inline double
pga2_scalar_product(pga2 a, pga2 b)
{
    const double *x = a.data, *y = b.data;
    return x[0]*y[0] + x[2]*y[2] + x[4]*y[4] - x[6]*y[6];
}

/* Write a = w + t*I + u, with I = e0*e1*e2 and I*I = 0.
 * The square u*u = q0 + q7*I is central. */
static inline double
pga2_quadratic_scalar(pga2 a)
{
    const double *x = a.data;
    return x[2]*x[2] + x[4]*x[4] - x[6]*x[6];
}

static inline double
pga2_quadratic_center(pga2 a)
{
    const double *x = a.data;
    return 2*(x[1]*x[6] - x[2]*x[5] + x[4]*x[3]);
}

static inline pga2
pga2_conjugate_product(pga2 a)
{
    double d0 = a.data[0]*a.data[0] - pga2_quadratic_scalar(a);
    double d7 = 2*a.data[0]*a.data[7] - pga2_quadratic_center(a);
    return pga2_make((const double[]){d0, 0, 0, 0, 0, 0, 0, d7});
}

static inline double
pga2_norm4(pga2 a)
{
    pga2 d = pga2_conjugate_product(a);
    return d.data[0]*d.data[0];
}

static inline double
pga2_norm(pga2 a)
{
    return sqrt(fabs(pga2_conjugate_product(a).data[0]));
}

static inline pga2
pga2_inv(pga2 a)
{
    pga2 d = pga2_conjugate_product(a);
    double norm = d.data[0]*d.data[0];
    pga2 dinv = pga2_make((const double[]){d.data[0], 0, 0, 0, 0, 0, 0, -d.data[7]});
    return vec8_div_scalar(pga2_mul(pga2_conjugate(a), dinv), norm);
}

static inline void
pga2_minpoly(pga2 a, double result[5])
{
    double w = a.data[0];
    double q0 = pga2_quadratic_scalar(a);
    result[0] = pga2_norm4(a);
    result[1] = -4*w*w*w + 4*w*q0;
    result[2] = 6*w*w - 2*q0;
    result[3] = -4*w;
    result[4] = 1;
}

static inline void
pga2_eivals(pga2 a, complex result[4])
{
    complex root = complex_sqrt_real(pga2_quadratic_scalar(a));
    result[0] = complex_make((const double[]){a.data[0] - root.data[0], -root.data[1]});
    result[1] = complex_make((const double[]){a.data[0] + root.data[0], root.data[1]});
    result[2] = result[0];
    result[3] = result[1];
}

static inline pga2
pga2_div(pga2 a, pga2 b)
{
    return pga2_mul(a, pga2_inv(b));
}

static inline pga2
pga2_real_slice(complex value)
{
    pga2 result = pga2_from_re(value.data[0]);
    result.data[6] = value.data[1];
    return result;
}

static inline pga2
pga2_center_slice(pga2 a, double a0, double a7, double b0, double b7)
{
    const double *x = a.data;
    return pga2_make((const double[]){a0, b0*x[1] - b7*x[6], b0*x[2], b0*x[3] + b7*x[4], b0*x[4], b0*x[5] - b7*x[2], b0*x[6], a7});
}

static inline int
pga2_is_central(pga2 a)
{
    const double *x = a.data;
    return x[1] == 0 && x[2] == 0 && x[3] == 0 &&
        x[4] == 0 && x[5] == 0 && x[6] == 0;
}

static inline pga2
pga2_apply_real_quadratic(pga2 a, double q,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    if (hermite_function_distance(complex_make((const double[]){a.data[0], 0}), f) < 0) {
        complex root = complex_sqrt_real(q);
        complex nodes[2] = {vec2_neg(root), root};
        double coefficients[4];
        if (hermite_interpolate(2, a.data[0], nodes,
                f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients)) {
            return pga2_center_slice(a, coefficients[0], 0, coefficients[1], 0);
        }
        return pga2_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
    }
    hermite_quadratic value = {0};
    if (hermite_quadratic_apply(a.data[0], q, 0,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return pga2_center_slice(a, value.even, 0, value.odd, 0);
    }
    return pga2_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

/* Distinct real channels can have a semisimple zero eigenvalue. Only
 * request a derivative when its projected nilpotent part is nonzero. */
static inline pga2
pga2_apply_real_channels(pga2 a, double q0, double q7,
    double (*f)(double), double (*deriv1)(double))
{
    double radius = sqrt(q0), shift = q7/(2*radius), t = a.data[7];
    double plus = f(a.data[0] + radius), minus = f(a.data[0] - radius);
    double odd = (plus/2 - minus/2)/radius;
    pga2 np = pga2_center_slice(a, 0, (t + shift)/2, 0, (t + shift)/(2*radius));
    pga2 nm = pga2_center_slice(a, 0, (t - shift)/2, 0, -(t - shift)/(2*radius));
    double dp = vec8_is_real(np) ? 0 : deriv1(a.data[0] + radius);
    double dm = vec8_is_real(nm) ? 0 : deriv1(a.data[0] - radius);
    if (isfinite(plus) && isfinite(minus) && isfinite(dp) && isfinite(dm)) {
        pga2 result = pga2_center_slice(a, plus/2 + minus/2, 0, odd, -q7*odd/(2*q0));
        return vec8_add(result, vec8_add(vec8_scale(np, dp), vec8_scale(nm, dm)));
    }
    return pga2_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline pga2
pga2_apply_dual_center(pga2 a, double q0, double q7,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    double t = a.data[7];
    if (q0 > 0) {
        double shift = q7/(2*sqrt(q0));
        if (t + shift == 0 || t - shift == 0) {
            return pga2_apply_real_channels(a, q0, q7, f, deriv1);
        }
    }
    if (hermite_quadratic_apply(a.data[0], q0, 1,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return pga2_center_slice(a, value.even,
            t*value.even_dw + q7*value.odd_dw/2, value.odd,
            t*value.odd_dw + q7*value.odd_dq);
    }
    if (q0 > 0) {
        return pga2_apply_real_channels(a, q0, q7, f, deriv1);
    }
    return pga2_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

/* The generic polynomial is (z*z - q0)^2 after subtracting the scalar.
 * Checked Hermite differences keep custom callbacks stable near coalescing
 * roots without adding that work to the elementary-function path. */
static inline pga2
pga2_apply_custom(pga2 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double w = a.data[0], coefficients[4];
    complex nodes[4];
    pga2 b = a;
    b.data[0] = 0;
    pga2_eivals(b, nodes);
    if (hermite_interpolate(4, w, nodes,
            f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients)) {
        pga2 b2 = pga2_square(b), b3 = pga2_mul(b2, b);
        pga2 result = vec8_add(vec8_scale(b, coefficients[1]),
            vec8_add(vec8_scale(b2, coefficients[2]), vec8_scale(b3, coefficients[3])));
        result.data[0] += coefficients[0];
        return result;
    }
    if (pga2_quadratic_scalar(a) > 0) {
        return pga2_apply_real_channels(a, pga2_quadratic_scalar(a),
            pga2_quadratic_center(a), f, deriv1);
    }
    return pga2_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline pga2
pga2_apply_central(pga2 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    (void)deriv2; (void)deriv3; (void)cf; (void)cderiv1;
    double value = f(a.data[0]), derivative = deriv1(a.data[0]);
    if (isfinite(value) && isfinite(derivative)) {
        return pga2_make((const double[]){value, 0, 0, 0, 0, 0, 0, a.data[7]*derivative});
    }
    return pga2_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline pga2
pga2_apply(pga2 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double q0 = pga2_quadratic_scalar(a), q7 = pga2_quadratic_center(a);
    if (vec8_is_real(a)) {
        return pga2_real_slice(cf(complex_make((const double[]){a.data[0], 0})));
    }
    if (a.data[7] == 0 && q7 == 0) {
        return pga2_apply_real_quadratic(a, q0, f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    if (pga2_is_central(a)) {
        return pga2_apply_central(a, f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    if (hermite_function_distance(complex_make((const double[]){a.data[0], 0}), f) < 0) {
        return pga2_apply_custom(a, f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    return pga2_apply_dual_center(a, q0, q7, f, deriv1, deriv2, deriv3, cf, cderiv1);
}

static inline pga2
pga2_exp(pga2 a)
{
    return pga2_apply(a, double_exp, double_exp_deriv1,
        double_exp_deriv2, double_exp_deriv3, complex_exp, complex_exp_deriv1);
}

static inline pga2
pga2_log(pga2 a)
{
    return pga2_apply(a, double_log, double_log_deriv1,
        double_log_deriv2, double_log_deriv3, complex_log, complex_log_deriv1);
}

static inline pga2
pga2_sqrt(pga2 a)
{
    return pga2_apply(a, double_sqrt, double_sqrt_deriv1,
        double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt, complex_sqrt_deriv1);
}

static inline pga2
pga2_cos(pga2 a)
{
    return pga2_apply(a, double_cos, double_cos_deriv1,
        double_cos_deriv2, double_cos_deriv3, complex_cos, complex_cos_deriv1);
}

static inline pga2
pga2_sin(pga2 a)
{
    return pga2_apply(a, double_sin, double_sin_deriv1,
        double_sin_deriv2, double_sin_deriv3, complex_sin, complex_sin_deriv1);
}

static inline pga2
pga2_tan(pga2 a)
{
    return pga2_apply(a, double_tan, double_tan_deriv1,
        double_tan_deriv2, double_tan_deriv3, complex_tan, complex_tan_deriv1);
}

static inline pga2
pga2_acos(pga2 a)
{
    return pga2_apply(a, double_acos, double_acos_deriv1,
        double_acos_deriv2, double_acos_deriv3, complex_acos, complex_acos_deriv1);
}

static inline pga2
pga2_asin(pga2 a)
{
    return pga2_apply(a, double_asin, double_asin_deriv1,
        double_asin_deriv2, double_asin_deriv3, complex_asin, complex_asin_deriv1);
}

static inline pga2
pga2_atan(pga2 a)
{
    return pga2_apply(a, double_atan, double_atan_deriv1,
        double_atan_deriv2, double_atan_deriv3, complex_atan, complex_atan_deriv1);
}

static inline pga2
pga2_cosh(pga2 a)
{
    return pga2_apply(a, double_cosh, double_cosh_deriv1,
        double_cosh_deriv2, double_cosh_deriv3, complex_cosh, complex_cosh_deriv1);
}

static inline pga2
pga2_sinh(pga2 a)
{
    return pga2_apply(a, double_sinh, double_sinh_deriv1,
        double_sinh_deriv2, double_sinh_deriv3, complex_sinh, complex_sinh_deriv1);
}

static inline pga2
pga2_tanh(pga2 a)
{
    return pga2_apply(a, double_tanh, double_tanh_deriv1,
        double_tanh_deriv2, double_tanh_deriv3, complex_tanh, complex_tanh_deriv1);
}

static inline pga2
pga2_acosh(pga2 a)
{
    return pga2_apply(a, double_acosh, double_acosh_deriv1,
        double_acosh_deriv2, double_acosh_deriv3, complex_acosh, complex_acosh_deriv1);
}

static inline pga2
pga2_asinh(pga2 a)
{
    return pga2_apply(a, double_asinh, double_asinh_deriv1,
        double_asinh_deriv2, double_asinh_deriv3, complex_asinh, complex_asinh_deriv1);
}

static inline pga2
pga2_atanh(pga2 a)
{
    return pga2_apply(a, double_atanh, double_atanh_deriv1,
        double_atanh_deriv2, double_atanh_deriv3, complex_atanh, complex_atanh_deriv1);
}

static inline pga2
pga2_expm1(pga2 a)
{
    return pga2_apply(a, double_expm1, double_expm1_deriv1,
        double_expm1_deriv2, double_expm1_deriv3, complex_expm1, complex_expm1_deriv1);
}

static inline pga2
pga2_log1p(pga2 a)
{
    return pga2_apply(a, double_log1p, double_log1p_deriv1,
        double_log1p_deriv2, double_log1p_deriv3, complex_log1p, complex_log1p_deriv1);
}

static inline pga2
pga2_radd(double a, pga2 b)
{
    return pga2_add(pga2_from_re(a), b);
}

static inline pga2
pga2_ladd(pga2 a, double b)
{
    return pga2_add(a, pga2_from_re(b));
}

static inline pga2
pga2_rsub(double a, pga2 b)
{
    return pga2_sub(pga2_from_re(a), b);
}

static inline pga2
pga2_lsub(pga2 a, double b)
{
    return pga2_sub(a, pga2_from_re(b));
}

static inline pga2
pga2_rmul(double a, pga2 b)
{
    return pga2_mul(pga2_from_re(a), b);
}

static inline pga2
pga2_lmul(pga2 a, double b)
{
    return pga2_mul(a, pga2_from_re(b));
}

static inline pga2
pga2_rdiv(double a, pga2 b)
{
    return pga2_div(pga2_from_re(a), b);
}

static inline pga2
pga2_ldiv(pga2 a, double b)
{
    return pga2_div(a, pga2_from_re(b));
}

static inline pga2
pga2_ipow(pga2 a, int64_t b)
{
    pga2 result = pga2_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return pga2_mul(a, a);
    }
    if (b == -1) {
        return pga2_inv(a);
    }
    if (b < 0) {
        a = pga2_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = pga2_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = pga2_mul(a, a);
        }
    }
    return result;
}

static inline pga2
pga2_rpow(pga2 a, pga2 b)
{
    return pga2_exp(pga2_mul(pga2_log(a), b));
}

static inline pga2
pga2_pow(pga2 a, pga2 b)
{
    return pga2_rpow(a, b);
}

static inline pga2
pga2_lpow(pga2 a, pga2 b)
{
    return pga2_exp(pga2_mul(b, pga2_log(a)));
}

static inline pga2
pga2_rpows(double a, pga2 b)
{
    return pga2_rpow(pga2_from_re(a), b);
}

static inline pga2
pga2_lpows(pga2 a, double b)
{
    return pga2_rpow(a, pga2_from_re(b));
}

static inline pga2
pga2_mix(pga2 a, pga2 b, double t)
{
    pga2 delta = pga2_mul(pga2_inv(a), b);
    return pga2_mul(a, pga2_lpows(delta, t));
}

static inline void
pga2_lmat(pga2 a, double result[ALGPY_PGA2_N*ALGPY_PGA2_N])
{
    pga2 basis = {{0}};
    pga2 product;
    int row, col;
    for (col = 0; col < ALGPY_PGA2_N; col++) {
        basis.data[col] = 1;
        product = pga2_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_PGA2_N; row++) {
            result[row * ALGPY_PGA2_N + col] = product.data[row];
        }
    }
}

static inline void
pga2_rmat(pga2 a, double result[ALGPY_PGA2_N*ALGPY_PGA2_N])
{
    pga2 basis = {{0}};
    pga2 product;
    int row, col;
    for (row = 0; row < ALGPY_PGA2_N; row++) {
        basis.data[row] = 1;
        product = pga2_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_PGA2_N; col++) {
            result[row * ALGPY_PGA2_N + col] = product.data[col];
        }
    }
}

/* Lattice API. The exterior and regressive products and the Poincare duals
 * do not depend on the metric, so they are those of exterior3, which shares
 * this coefficient order: join with identity bottom = 1, meet with identity
 * top = the pseudoscalar, rdual(e_a) join e_a = top = e_a join ldual(e_a). */
static inline pga2
pga2_bottom(void)
{
    return exterior3_bottom();
}

static inline pga2
pga2_top(void)
{
    return exterior3_top();
}

static inline pga2
pga2_join(pga2 a, pga2 b)
{
    return exterior3_join(a, b);
}

static inline pga2
pga2_meet(pga2 a, pga2 b)
{
    return exterior3_meet(a, b);
}

static inline pga2
pga2_ldual(pga2 a)
{
    return exterior3_ldual(a);
}

static inline pga2
pga2_rdual(pga2 a)
{
    return exterior3_rdual(a);
}

static inline pga2
pga2_dual(pga2 a)
{
    return exterior3_dual(a);
}

/*
 * Minimal polynomial over the centre span{1, w}, w = e012 with w*w = 0,
 * which is the dual1 algebra. Coefficients are in ascending order, p[i]
 * standing for p[i][0] + p[i][1]*w: the central quadratic
 * a*a - (a + conj(a))*a + a*conj(a), with a + conj(a) = 2*(x + x012*w).
 */
static inline void
pga2_center_minpoly(pga2 a, dual1 p[3])
{
    pga2 delta = pga2_conjugate_product(a);
    p[0] = dual1_make((const double[]){delta.data[0], delta.data[7]});
    p[1] = dual1_make((const double[]){-2*a.data[0], -2*a.data[7]});
    p[2] = dual1_make((const double[]){1, 0});
}

#endif /* ALGPY_PGA2_H */
