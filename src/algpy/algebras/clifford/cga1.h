/* Conformal geometric algebra in 1D.
 * Generator squares: -1, 1, 1.
 * 3 generators; 8 real coefficients. */
#ifndef ALGPY_CGA1_H
#define ALGPY_CGA1_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../polynomials/hermite.h"
#include "../../vectorspaces/vec8.h"
#include "../../lattices/clifford/exterior3.h"
#include "../multicomplex/splitcomplex1.h"

#define ALGPY_CGA1_N 8

typedef vec8 cga1;

static inline cga1
cga1_make(const double coeffs[ALGPY_CGA1_N])
{
    cga1 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline cga1
cga1_grade0(cga1 a)
{
    const double *x = a.data;
    return cga1_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0});
}

static inline cga1
cga1_grade1(cga1 a)
{
    const double *x = a.data;
    return cga1_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0});
}

static inline cga1
cga1_grade2(cga1 a)
{
    const double *x = a.data;
    return cga1_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0});
}

static inline cga1
cga1_grade3(cga1 a)
{
    const double *x = a.data;
    return cga1_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7]});
}

static inline cga1
cga1_from_re(double a)
{
    return cga1_make((const double[ALGPY_CGA1_N]){a});
}

static inline double
cga1_to_re(cga1 a)
{
    return a.data[0];
}

static inline cga1
cga1_zero(void)
{
    return vec8_zero();
}

static inline cga1
cga1_one(void)
{
    return cga1_from_re(1);
}

static inline cga1
cga1_pos(cga1 a)
{
    return vec8_pos(a);
}

static inline cga1
cga1_neg(cga1 a)
{
    return vec8_neg(a);
}

static inline cga1
cga1_add(cga1 a, cga1 b)
{
    return vec8_add(a, b);
}

static inline cga1
cga1_sub(cga1 a, cga1 b)
{
    return vec8_sub(a, b);
}

static inline int
cga1_equal(cga1 a, cga1 b)
{
    return vec8_all_eq(a, b);
}

static inline int
cga1_not_equal(cga1 a, cga1 b)
{
    return !cga1_equal(a, b);
}

static inline int
cga1_is_real(cga1 a)
{
    return vec8_is_real(a);
}

static inline cga1
cga1_re(cga1 a)
{
    return cga1_from_re(a.data[0]);
}

static inline cga1
cga1_im(cga1 a)
{
    return cga1_sub(a, cga1_re(a));
}

static inline cga1
cga1_mul(cga1 a, cga1 b)
{
    const double *x = a.data, *y = b.data;
    return cga1_make((const double[]){
        x[0]*y[0] - x[1]*y[1] + x[2]*y[2] + x[3]*y[3] + x[4]*y[4] + x[5]*y[5] - x[6]*y[6] + x[7]*y[7],
        x[0]*y[1] + x[1]*y[0] - x[2]*y[3] + x[3]*y[2] - x[4]*y[5] + x[5]*y[4] - x[6]*y[7] - x[7]*y[6],
        x[0]*y[2] - x[1]*y[3] + x[2]*y[0] + x[3]*y[1] - x[4]*y[6] - x[5]*y[7] + x[6]*y[4] - x[7]*y[5],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0] + x[4]*y[7] - x[5]*y[6] + x[6]*y[5] + x[7]*y[4],
        x[0]*y[4] - x[1]*y[5] + x[2]*y[6] + x[3]*y[7] + x[4]*y[0] + x[5]*y[1] - x[6]*y[2] + x[7]*y[3],
        x[0]*y[5] + x[1]*y[4] - x[2]*y[7] + x[3]*y[6] - x[4]*y[1] + x[5]*y[0] - x[6]*y[3] - x[7]*y[2],
        x[0]*y[6] - x[1]*y[7] + x[2]*y[4] + x[3]*y[5] - x[4]*y[2] - x[5]*y[3] + x[6]*y[0] - x[7]*y[1],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0]});
}

static inline cga1
cga1_square(cga1 a)
{
    const double *x = a.data;
    return cga1_make((const double[]){
        x[0]*x[0] - x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5] - x[6]*x[6] + x[7]*x[7],
        2*x[0]*x[1] - 2*x[6]*x[7],
        2*x[0]*x[2] - 2*x[5]*x[7],
        2*x[0]*x[3] + 2*x[4]*x[7],
        2*x[0]*x[4] + 2*x[3]*x[7],
        2*x[0]*x[5] - 2*x[2]*x[7],
        2*x[0]*x[6] - 2*x[1]*x[7],
        2*x[0]*x[7] + 2*x[1]*x[6] - 2*x[2]*x[5] + 2*x[3]*x[4]});
}

static inline cga1
cga1_conjugate(cga1 a)
{
    const double *x = a.data;
    return cga1_make((const double[]){x[0], -x[1], -x[2], -x[3], -x[4], -x[5], -x[6], x[7]});
}

static inline void
cga1_basis(cga1 result[8])
{
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cga1_gens(cga1 result[3])
{
    result[0] = cga1_make((const double[ALGPY_CGA1_N]){[1] = 1});
    result[1] = cga1_make((const double[ALGPY_CGA1_N]){[2] = 1});
    result[2] = cga1_make((const double[ALGPY_CGA1_N]){[4] = 1});
}

static inline void
cga1_coeffs(cga1 result[64])
{
    cga1 basis[8];
    int i, j;
    cga1_basis(basis);
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            result[i*8+j] = cga1_mul(basis[i], basis[j]);
        }
    }
}

static inline double
cga1_scalar_product(cga1 a, cga1 b)
{
    const double *x = a.data, *y = b.data;
    return x[0]*y[0] - x[1]*y[1] + x[2]*y[2] + x[3]*y[3] + x[4]*y[4] + x[5]*y[5] - x[6]*y[6] + x[7]*y[7];
}

/* Write a = w + t*I + u, with I = e0*e1*e2 and I*I = 1.
 * The square u*u = q0 + q7*I splits into q0 +/- q7. */
static inline double
cga1_quadratic_scalar(cga1 a)
{
    const double *x = a.data;
    return -x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5] - x[6]*x[6];
}

static inline double
cga1_quadratic_center(cga1 a)
{
    const double *x = a.data;
    return 2*(x[1]*x[6] - x[2]*x[5] + x[4]*x[3]);
}

static inline cga1
cga1_conjugate_product(cga1 a)
{
    double d0 = a.data[0]*a.data[0] + a.data[7]*a.data[7] - cga1_quadratic_scalar(a);
    double d7 = 2*a.data[0]*a.data[7] - cga1_quadratic_center(a);
    return cga1_make((const double[]){d0, 0, 0, 0, 0, 0, 0, d7});
}

static inline double
cga1_norm4(cga1 a)
{
    cga1 d = cga1_conjugate_product(a);
    return d.data[0]*d.data[0] - d.data[7]*d.data[7];
}

static inline complex
cga1_norm(cga1 a)
{
    return complex_root_real(cga1_norm4(a), 4);
}

static inline cga1
cga1_inv(cga1 a)
{
    cga1 d = cga1_conjugate_product(a);
    double norm = d.data[0]*d.data[0] - d.data[7]*d.data[7];
    cga1 dinv = cga1_make((const double[]){d.data[0], 0, 0, 0, 0, 0, 0, -d.data[7]});
    return vec8_div_scalar(cga1_mul(cga1_conjugate(a), dinv), norm);
}

static inline void
cga1_minpoly(cga1 a, double result[5])
{
    double w = a.data[0], t = a.data[7];
    double q0 = cga1_quadratic_scalar(a), q7 = cga1_quadratic_center(a);
    double s2 = q0 + t*t;
    result[0] = cga1_norm4(a);
    result[1] = -4*w*w*w + 4*w*s2 - 4*t*q7;
    result[2] = 6*w*w - 2*s2;
    result[3] = -4*w;
    result[4] = 1;
}

static inline void
cga1_eivals(cga1 a, complex result[4])
{
    double q0 = cga1_quadratic_scalar(a), q7 = cga1_quadratic_center(a);
    complex minus = complex_sqrt_real(q0 - q7);
    complex plus = complex_sqrt_real(q0 + q7);
    double wm = a.data[0] - a.data[7], wp = a.data[0] + a.data[7];
    result[0] = complex_make((const double[]){wm - minus.data[0], -minus.data[1]});
    result[1] = complex_make((const double[]){wm + minus.data[0], minus.data[1]});
    result[2] = complex_make((const double[]){wp - plus.data[0], -plus.data[1]});
    result[3] = complex_make((const double[]){wp + plus.data[0], plus.data[1]});
}

static inline cga1
cga1_div(cga1 a, cga1 b)
{
    return cga1_mul(a, cga1_inv(b));
}

static inline cga1
cga1_real_slice(complex value)
{
    cga1 result = cga1_from_re(value.data[0]);
    result.data[1] = value.data[1];
    return result;
}

static inline cga1
cga1_center_slice(cga1 a, double a0, double a7, double b0, double b7)
{
    const double *x = a.data;
    return cga1_make((const double[]){a0, b0*x[1] - b7*x[6], b0*x[2] - b7*x[5], b0*x[3] + b7*x[4], b0*x[4] + b7*x[3], b0*x[5] - b7*x[2], b0*x[6] - b7*x[1], a7});
}

static inline int
cga1_is_central(cga1 a)
{
    const double *x = a.data;
    return x[1] == 0 && x[2] == 0 && x[3] == 0 &&
        x[4] == 0 && x[5] == 0 && x[6] == 0;
}

/* Real or conjugate quadratic channels; custom callbacks use checked
 * Hermite interpolation when the roots coalesce. */
static inline int
cga1_quadratic(double w, double q,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), hermite_quadratic *result)
{
    /* The interpolation engine checks derivative quadrature for custom
     * callbacks, avoiding cancellation when a channel has repeated roots. */
    if (hermite_function_distance(complex_make((const double[]){w, 0}), f) < 0) {
        complex root = complex_sqrt_real(q);
        complex nodes[2] = {vec2_neg(root), root};
        double coefficients[4];
        if (!hermite_interpolate(2, w, nodes,
                f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients)) {
            return 0;
        }
        result->even = coefficients[0];
        result->odd = coefficients[1];
        return 1;
    }
    return hermite_quadratic_apply(w, q, 0,
        f, deriv1, deriv2, deriv3, cf, cderiv1, result);
}

static inline cga1
cga1_apply_real_quadratic(cga1 a, double q,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    if (cga1_quadratic(a.data[0], q,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return cga1_center_slice(a, value.even, 0, value.odd, 0);
    }
    return cga1_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

typedef struct {
    double even_re, even_i, odd_re, odd_i;
} cga1_center_coefficients;

static inline int
cga1_center_channels_are_close(double w, double t, double q0, double q7,
    double (*f)(double))
{
    double radius, distance;
    radius = sqrt(fabs(q0));
    if (q0 >= 0) {
        distance = fmin(hermite_function_distance(complex_make((const double[]){w + radius, 0}), f),
            hermite_function_distance(complex_make((const double[]){w - radius, 0}), f));
    } else {
        distance = hermite_function_distance(complex_make((const double[]){w, radius}), f);
    }
    return distance > 0 && (fabs(t)/distance +
        fabs(q7)/distance/distance) < 0.001 && isfinite(f(w));
}

static inline int
cga1_near_center_coefficients(double w, double t, double q0, double q7,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), cga1_center_coefficients *result)
{
    hermite_quadratic base = {0};
    static const double odd[4] = {586./945, 167./420, -2./105, 1./756};
    static const double even[3] = {179./576, -13./360, 11./2880};
    int i;
    if (!hermite_quadratic_apply(w, q0, 1,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &base)) {
        return 0;
    }
    result->even_re = base.even;
    result->even_i = odd[0]*(t*base.even_dw + q7*base.odd_dw/2);
    result->odd_re = base.odd;
    result->odd_i = odd[0]*(t*base.odd_dw + q7*base.odd_dq);
    for (i = 1; i <= 3; i++) {
        double offset = i;
        hermite_quadratic plus = {0}, minus = {0};
        double ap, am, bp, bm;
        if (!hermite_quadratic_apply(w + t*offset, q0 + q7*offset, 1,
                f, deriv1, deriv2, deriv3, cf, cderiv1, &plus) ||
            !hermite_quadratic_apply(w - t*offset, q0 - q7*offset, 1,
                f, deriv1, deriv2, deriv3, cf, cderiv1, &minus)) {
            return 0;
        }
        ap = t*plus.even_dw + q7*plus.odd_dw/2;
        am = t*minus.even_dw + q7*minus.odd_dw/2;
        bp = t*plus.odd_dw + q7*plus.odd_dq;
        bm = t*minus.odd_dw + q7*minus.odd_dq;
        result->even_re += even[i - 1]*(ap - am);
        result->even_i += odd[i]*(ap/2 + am/2);
        result->odd_re += even[i - 1]*(bp - bm);
        result->odd_i += odd[i]*(bp/2 + bm/2);
    }
    return 1;
}

static inline int
cga1_split_channel(cga1 u, double w, double q,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), hermite_quadratic *result)
{
    if (vec8_is_real(u)) {
        result->even = f(w);
        result->odd = 0;
        return isfinite(result->even);
    }
    return cga1_quadratic(w, q,
        f, deriv1, deriv2, deriv3, cf, cderiv1, result);
}

static inline cga1
cga1_apply_split_center(cga1 a, double q0, double q7,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic plus = {0}, minus = {0};
    cga1 up = cga1_center_slice(a, 0, 0, 1, 1);
    cga1 um = cga1_center_slice(a, 0, 0, 1, -1);
    int plus_ok = cga1_split_channel(up, a.data[0] + a.data[7], q0 + q7,
        f, deriv1, deriv2, deriv3, cf, cderiv1, &plus);
    int minus_ok = cga1_split_channel(um, a.data[0] - a.data[7], q0 - q7,
        f, deriv1, deriv2, deriv3, cf, cderiv1, &minus);
    if (plus_ok && minus_ok) {
        return cga1_center_slice(a,
            plus.even/2 + minus.even/2, (plus.even/2 - minus.even/2),
            plus.odd/2 + minus.odd/2, (plus.odd/2 - minus.odd/2));
    }
    return cga1_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline cga1
cga1_apply_central(cga1 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    if (cga1_quadratic(a.data[0], a.data[7]*a.data[7],
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return cga1_make((const double[]){value.even, 0, 0, 0, 0, 0, 0, a.data[7]*value.odd});
    }
    return cga1_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline cga1
cga1_apply(cga1 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double q0 = cga1_quadratic_scalar(a), q7 = cga1_quadratic_center(a);
    if (vec8_is_real(a)) {
        return cga1_real_slice(cf(complex_make((const double[]){a.data[0], 0})));
    }
    if (a.data[7] == 0 && q7 == 0) {
        return cga1_apply_real_quadratic(a, q0, f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    if (cga1_is_central(a)) {
        return cga1_apply_central(a, f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    if (cga1_center_channels_are_close(a.data[0], a.data[7], q0, q7, f)) {
        cga1_center_coefficients value;
        if (cga1_near_center_coefficients(a.data[0], a.data[7], q0, q7,
                f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
            return cga1_center_slice(a, value.even_re, value.even_i,
                value.odd_re, value.odd_i);
        }
    }
    return cga1_apply_split_center(a, q0, q7, f, deriv1, deriv2, deriv3, cf, cderiv1);
}

static inline cga1
cga1_exp(cga1 a)
{
    return cga1_apply(a, double_exp, double_exp_deriv1,
        double_exp_deriv2, double_exp_deriv3, complex_exp, complex_exp_deriv1);
}

static inline cga1
cga1_log(cga1 a)
{
    return cga1_apply(a, double_log, double_log_deriv1,
        double_log_deriv2, double_log_deriv3, complex_log, complex_log_deriv1);
}

static inline cga1
cga1_sqrt(cga1 a)
{
    return cga1_apply(a, double_sqrt, double_sqrt_deriv1,
        double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt, complex_sqrt_deriv1);
}

static inline cga1
cga1_cos(cga1 a)
{
    return cga1_apply(a, double_cos, double_cos_deriv1,
        double_cos_deriv2, double_cos_deriv3, complex_cos, complex_cos_deriv1);
}

static inline cga1
cga1_sin(cga1 a)
{
    return cga1_apply(a, double_sin, double_sin_deriv1,
        double_sin_deriv2, double_sin_deriv3, complex_sin, complex_sin_deriv1);
}

static inline cga1
cga1_tan(cga1 a)
{
    return cga1_apply(a, double_tan, double_tan_deriv1,
        double_tan_deriv2, double_tan_deriv3, complex_tan, complex_tan_deriv1);
}

static inline cga1
cga1_acos(cga1 a)
{
    return cga1_apply(a, double_acos, double_acos_deriv1,
        double_acos_deriv2, double_acos_deriv3, complex_acos, complex_acos_deriv1);
}

static inline cga1
cga1_asin(cga1 a)
{
    return cga1_apply(a, double_asin, double_asin_deriv1,
        double_asin_deriv2, double_asin_deriv3, complex_asin, complex_asin_deriv1);
}

static inline cga1
cga1_atan(cga1 a)
{
    return cga1_apply(a, double_atan, double_atan_deriv1,
        double_atan_deriv2, double_atan_deriv3, complex_atan, complex_atan_deriv1);
}

static inline cga1
cga1_cosh(cga1 a)
{
    return cga1_apply(a, double_cosh, double_cosh_deriv1,
        double_cosh_deriv2, double_cosh_deriv3, complex_cosh, complex_cosh_deriv1);
}

static inline cga1
cga1_sinh(cga1 a)
{
    return cga1_apply(a, double_sinh, double_sinh_deriv1,
        double_sinh_deriv2, double_sinh_deriv3, complex_sinh, complex_sinh_deriv1);
}

static inline cga1
cga1_tanh(cga1 a)
{
    return cga1_apply(a, double_tanh, double_tanh_deriv1,
        double_tanh_deriv2, double_tanh_deriv3, complex_tanh, complex_tanh_deriv1);
}

static inline cga1
cga1_acosh(cga1 a)
{
    return cga1_apply(a, double_acosh, double_acosh_deriv1,
        double_acosh_deriv2, double_acosh_deriv3, complex_acosh, complex_acosh_deriv1);
}

static inline cga1
cga1_asinh(cga1 a)
{
    return cga1_apply(a, double_asinh, double_asinh_deriv1,
        double_asinh_deriv2, double_asinh_deriv3, complex_asinh, complex_asinh_deriv1);
}

static inline cga1
cga1_atanh(cga1 a)
{
    return cga1_apply(a, double_atanh, double_atanh_deriv1,
        double_atanh_deriv2, double_atanh_deriv3, complex_atanh, complex_atanh_deriv1);
}

static inline cga1
cga1_expm1(cga1 a)
{
    return cga1_apply(a, double_expm1, double_expm1_deriv1,
        double_expm1_deriv2, double_expm1_deriv3, complex_expm1, complex_expm1_deriv1);
}

static inline cga1
cga1_log1p(cga1 a)
{
    return cga1_apply(a, double_log1p, double_log1p_deriv1,
        double_log1p_deriv2, double_log1p_deriv3, complex_log1p, complex_log1p_deriv1);
}

static inline cga1
cga1_radd(double a, cga1 b)
{
    return cga1_add(cga1_from_re(a), b);
}

static inline cga1
cga1_ladd(cga1 a, double b)
{
    return cga1_add(a, cga1_from_re(b));
}

static inline cga1
cga1_rsub(double a, cga1 b)
{
    return cga1_sub(cga1_from_re(a), b);
}

static inline cga1
cga1_lsub(cga1 a, double b)
{
    return cga1_sub(a, cga1_from_re(b));
}

static inline cga1
cga1_rmul(double a, cga1 b)
{
    return cga1_mul(cga1_from_re(a), b);
}

static inline cga1
cga1_lmul(cga1 a, double b)
{
    return cga1_mul(a, cga1_from_re(b));
}

static inline cga1
cga1_rdiv(double a, cga1 b)
{
    return cga1_div(cga1_from_re(a), b);
}

static inline cga1
cga1_ldiv(cga1 a, double b)
{
    return cga1_div(a, cga1_from_re(b));
}

static inline cga1
cga1_ipow(cga1 a, int64_t b)
{
    cga1 result = cga1_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cga1_mul(a, a);
    }
    if (b == -1) {
        return cga1_inv(a);
    }
    if (b < 0) {
        a = cga1_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = cga1_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = cga1_mul(a, a);
        }
    }
    return result;
}

static inline cga1
cga1_rpow(cga1 a, cga1 b)
{
    return cga1_exp(cga1_mul(cga1_log(a), b));
}

static inline cga1
cga1_pow(cga1 a, cga1 b)
{
    return cga1_rpow(a, b);
}

static inline cga1
cga1_lpow(cga1 a, cga1 b)
{
    return cga1_exp(cga1_mul(b, cga1_log(a)));
}

static inline cga1
cga1_rpows(double a, cga1 b)
{
    return cga1_rpow(cga1_from_re(a), b);
}

static inline cga1
cga1_lpows(cga1 a, double b)
{
    return cga1_rpow(a, cga1_from_re(b));
}

static inline cga1
cga1_mix(cga1 a, cga1 b, double t)
{
    cga1 delta = cga1_mul(cga1_inv(a), b);
    return cga1_mul(a, cga1_lpows(delta, t));
}

static inline void
cga1_lmat(cga1 a, double result[ALGPY_CGA1_N*ALGPY_CGA1_N])
{
    cga1 basis = {{0}};
    cga1 product;
    int row, col;
    for (col = 0; col < ALGPY_CGA1_N; col++) {
        basis.data[col] = 1;
        product = cga1_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CGA1_N; row++) {
            result[row * ALGPY_CGA1_N + col] = product.data[row];
        }
    }
}

static inline void
cga1_rmat(cga1 a, double result[ALGPY_CGA1_N*ALGPY_CGA1_N])
{
    cga1 basis = {{0}};
    cga1 product;
    int row, col;
    for (row = 0; row < ALGPY_CGA1_N; row++) {
        basis.data[row] = 1;
        product = cga1_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CGA1_N; col++) {
            result[row * ALGPY_CGA1_N + col] = product.data[col];
        }
    }
}

/* Lattice API. The exterior and regressive products and the Poincare duals
 * do not depend on the metric, so they are those of exterior3, which shares
 * this coefficient order: join with identity bottom = 1, meet with identity
 * top = the pseudoscalar, rdual(e_a) join e_a = top = e_a join ldual(e_a). */
static inline cga1
cga1_bottom(void)
{
    return exterior3_bottom();
}

static inline cga1
cga1_top(void)
{
    return exterior3_top();
}

static inline cga1
cga1_join(cga1 a, cga1 b)
{
    return exterior3_join(a, b);
}

static inline cga1
cga1_meet(cga1 a, cga1 b)
{
    return exterior3_meet(a, b);
}

static inline cga1
cga1_ldual(cga1 a)
{
    return exterior3_ldual(a);
}

static inline cga1
cga1_rdual(cga1 a)
{
    return exterior3_rdual(a);
}

static inline cga1
cga1_dual(cga1 a)
{
    return exterior3_dual(a);
}

/*
 * Minimal polynomial over the centre span{1, w}, w = e012 with w*w = 1,
 * which is the splitcomplex1 algebra. Coefficients are in ascending order, p[i]
 * standing for p[i][0] + p[i][1]*w: the central quadratic
 * a*a - (a + conj(a))*a + a*conj(a), with a + conj(a) = 2*(x + x012*w).
 */
static inline void
cga1_center_minpoly(cga1 a, splitcomplex1 p[3])
{
    cga1 delta = cga1_conjugate_product(a);
    p[0] = splitcomplex1_make((const double[]){delta.data[0], delta.data[7]});
    p[1] = splitcomplex1_make((const double[]){-2*a.data[0], -2*a.data[7]});
    p[2] = splitcomplex1_make((const double[]){1, 0});
}

#endif /* ALGPY_CGA1_H */
