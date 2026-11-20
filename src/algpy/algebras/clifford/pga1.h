/* Projective geometric algebra in 1D.
 * Generator squares: 0, 1.
 * 2 generators; 4 real coefficients. */
#ifndef ALGPY_PGA1_H
#define ALGPY_PGA1_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../polynomials/hermite.h"
#include "../../vectorspaces/vec4.h"
#include "../../lattices/clifford/exterior2.h"

#define ALGPY_PGA1_N 4

typedef vec4 pga1;

static inline pga1
pga1_make(const double coeffs[ALGPY_PGA1_N])
{
    pga1 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline pga1
pga1_grade0(pga1 a)
{
    const double *x = a.data;
    return pga1_make((const double[]){x[0], 0, 0, 0});
}

static inline pga1
pga1_grade1(pga1 a)
{
    const double *x = a.data;
    return pga1_make((const double[]){0, x[1], x[2], 0});
}

static inline pga1
pga1_grade2(pga1 a)
{
    const double *x = a.data;
    return pga1_make((const double[]){0, 0, 0, x[3]});
}

static inline pga1
pga1_from_re(double a)
{
    return pga1_make((const double[ALGPY_PGA1_N]){a});
}

static inline double
pga1_to_re(pga1 a)
{
    return a.data[0];
}

static inline pga1
pga1_zero(void)
{
    return vec4_zero();
}

static inline pga1
pga1_one(void)
{
    return pga1_from_re(1);
}

static inline pga1
pga1_pos(pga1 a)
{
    return vec4_pos(a);
}

static inline pga1
pga1_neg(pga1 a)
{
    return vec4_neg(a);
}

static inline pga1
pga1_add(pga1 a, pga1 b)
{
    return vec4_add(a, b);
}

static inline pga1
pga1_sub(pga1 a, pga1 b)
{
    return vec4_sub(a, b);
}

static inline int
pga1_equal(pga1 a, pga1 b)
{
    return vec4_all_eq(a, b);
}

static inline int
pga1_not_equal(pga1 a, pga1 b)
{
    return !pga1_equal(a, b);
}

static inline int
pga1_is_real(pga1 a)
{
    return vec4_is_real(a);
}

static inline pga1
pga1_re(pga1 a)
{
    return pga1_from_re(a.data[0]);
}

static inline pga1
pga1_im(pga1 a)
{
    return pga1_sub(a, pga1_re(a));
}

static inline pga1
pga1_mul(pga1 a, pga1 b)
{
    const double *x = a.data, *y = b.data;
    return pga1_make((const double[]){
        x[0]*y[0] + x[2]*y[2],
        x[0]*y[1] + x[1]*y[0] - x[2]*y[3] + x[3]*y[2],
        x[0]*y[2] + x[2]*y[0],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0]});
}

static inline pga1
pga1_square(pga1 a)
{
    const double *x = a.data;
    return pga1_make((const double[]){
        x[0]*x[0] + x[2]*x[2],
        2*x[0]*x[1],
        2*x[0]*x[2],
        2*x[0]*x[3]});
}

static inline pga1
pga1_conjugate(pga1 a)
{
    const double *x = a.data;
    return pga1_make((const double[]){x[0], -x[1], -x[2], -x[3]});
}

static inline void
pga1_basis(pga1 result[4])
{
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
pga1_gens(pga1 result[2])
{
    result[0] = pga1_make((const double[ALGPY_PGA1_N]){[1] = 1});
    result[1] = pga1_make((const double[ALGPY_PGA1_N]){[2] = 1});
}

static inline void
pga1_coeffs(pga1 result[16])
{
    pga1 basis[4];
    int i, j;
    pga1_basis(basis);
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            result[i*4+j] = pga1_mul(basis[i], basis[j]);
        }
    }
}

static inline double
pga1_scalar_product(pga1 a, pga1 b)
{
    const double *x = a.data, *y = b.data;
    return x[0]*y[0] + x[2]*y[2];
}

static inline double
pga1_norm2(pga1 a)
{
    return a.data[0]*a.data[0] - a.data[2]*a.data[2];
}

static inline complex
pga1_norm(pga1 a)
{
    return complex_sqrt_real(pga1_norm2(a));
}

static inline pga1
pga1_inv(pga1 a)
{
    return vec4_div_scalar(pga1_conjugate(a), pga1_norm2(a));
}

static inline void
pga1_minpoly(pga1 a, double result[3])
{
    result[0] = pga1_norm2(a);
    result[1] = -2*a.data[0];
    result[2] = 1;
}

static inline void
pga1_eivals(pga1 a, complex result[2])
{
    double radius = fabs(a.data[2]);
    result[0] = complex_make((const double[]){a.data[0] - radius, 0});
    result[1] = complex_make((const double[]){a.data[0] + radius, 0});
}

static inline pga1
pga1_div(pga1 a, pga1 b)
{
    return pga1_mul(a, pga1_inv(b));
}

static inline pga1
pga1_real_slice(complex value)
{
    pga1 result = pga1_from_re(value.data[0]);
    if (value.data[1] != 0) {
        return pga1_make((const double[]){NAN, NAN, NAN, NAN});
    }
    return result;
}

/* With a = w + v, v*v = a[2]*a[2]. The null coefficients do not
 * change the two roots, but are retained by the linear interpolant. */
static inline pga1
pga1_apply(pga1 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    if (vec4_is_real(a)) {
        return pga1_real_slice(cf(complex_make((const double[]){a.data[0], 0})));
    }
    if (hermite_function_distance(complex_make((const double[]){a.data[0], 0}), f) < 0) {
        double coefficients[4];
        complex nodes[2] = {complex_make((const double[]){-fabs(a.data[2]), 0}), complex_make((const double[]){fabs(a.data[2]), 0})};
        if (hermite_interpolate(2, a.data[0], nodes,
                f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients)) {
            return pga1_make((const double[]){coefficients[0], coefficients[1]*a.data[1],
                coefficients[1]*a.data[2], coefficients[1]*a.data[3]});
        }
        return pga1_make((const double[]){NAN, NAN, NAN, NAN});
    }
    if (hermite_quadratic_apply(a.data[0], a.data[2]*a.data[2], 0,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return pga1_make((const double[]){value.even, value.odd*a.data[1], value.odd*a.data[2], value.odd*a.data[3]});
    }
    return pga1_make((const double[]){NAN, NAN, NAN, NAN});
}

static inline pga1
pga1_exp(pga1 a)
{
    return pga1_apply(a, double_exp, double_exp_deriv1,
        double_exp_deriv2, double_exp_deriv3, complex_exp, complex_exp_deriv1);
}

static inline pga1
pga1_log(pga1 a)
{
    return pga1_apply(a, double_log, double_log_deriv1,
        double_log_deriv2, double_log_deriv3, complex_log, complex_log_deriv1);
}

static inline pga1
pga1_sqrt(pga1 a)
{
    return pga1_apply(a, double_sqrt, double_sqrt_deriv1,
        double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt, complex_sqrt_deriv1);
}

static inline pga1
pga1_cos(pga1 a)
{
    return pga1_apply(a, double_cos, double_cos_deriv1,
        double_cos_deriv2, double_cos_deriv3, complex_cos, complex_cos_deriv1);
}

static inline pga1
pga1_sin(pga1 a)
{
    return pga1_apply(a, double_sin, double_sin_deriv1,
        double_sin_deriv2, double_sin_deriv3, complex_sin, complex_sin_deriv1);
}

static inline pga1
pga1_tan(pga1 a)
{
    return pga1_apply(a, double_tan, double_tan_deriv1,
        double_tan_deriv2, double_tan_deriv3, complex_tan, complex_tan_deriv1);
}

static inline pga1
pga1_acos(pga1 a)
{
    return pga1_apply(a, double_acos, double_acos_deriv1,
        double_acos_deriv2, double_acos_deriv3, complex_acos, complex_acos_deriv1);
}

static inline pga1
pga1_asin(pga1 a)
{
    return pga1_apply(a, double_asin, double_asin_deriv1,
        double_asin_deriv2, double_asin_deriv3, complex_asin, complex_asin_deriv1);
}

static inline pga1
pga1_atan(pga1 a)
{
    return pga1_apply(a, double_atan, double_atan_deriv1,
        double_atan_deriv2, double_atan_deriv3, complex_atan, complex_atan_deriv1);
}

static inline pga1
pga1_cosh(pga1 a)
{
    return pga1_apply(a, double_cosh, double_cosh_deriv1,
        double_cosh_deriv2, double_cosh_deriv3, complex_cosh, complex_cosh_deriv1);
}

static inline pga1
pga1_sinh(pga1 a)
{
    return pga1_apply(a, double_sinh, double_sinh_deriv1,
        double_sinh_deriv2, double_sinh_deriv3, complex_sinh, complex_sinh_deriv1);
}

static inline pga1
pga1_tanh(pga1 a)
{
    return pga1_apply(a, double_tanh, double_tanh_deriv1,
        double_tanh_deriv2, double_tanh_deriv3, complex_tanh, complex_tanh_deriv1);
}

static inline pga1
pga1_acosh(pga1 a)
{
    return pga1_apply(a, double_acosh, double_acosh_deriv1,
        double_acosh_deriv2, double_acosh_deriv3, complex_acosh, complex_acosh_deriv1);
}

static inline pga1
pga1_asinh(pga1 a)
{
    return pga1_apply(a, double_asinh, double_asinh_deriv1,
        double_asinh_deriv2, double_asinh_deriv3, complex_asinh, complex_asinh_deriv1);
}

static inline pga1
pga1_atanh(pga1 a)
{
    return pga1_apply(a, double_atanh, double_atanh_deriv1,
        double_atanh_deriv2, double_atanh_deriv3, complex_atanh, complex_atanh_deriv1);
}

static inline pga1
pga1_expm1(pga1 a)
{
    return pga1_apply(a, double_expm1, double_expm1_deriv1,
        double_expm1_deriv2, double_expm1_deriv3, complex_expm1, complex_expm1_deriv1);
}

static inline pga1
pga1_log1p(pga1 a)
{
    return pga1_apply(a, double_log1p, double_log1p_deriv1,
        double_log1p_deriv2, double_log1p_deriv3, complex_log1p, complex_log1p_deriv1);
}

static inline pga1
pga1_radd(double a, pga1 b)
{
    return pga1_add(pga1_from_re(a), b);
}

static inline pga1
pga1_ladd(pga1 a, double b)
{
    return pga1_add(a, pga1_from_re(b));
}

static inline pga1
pga1_rsub(double a, pga1 b)
{
    return pga1_sub(pga1_from_re(a), b);
}

static inline pga1
pga1_lsub(pga1 a, double b)
{
    return pga1_sub(a, pga1_from_re(b));
}

static inline pga1
pga1_rmul(double a, pga1 b)
{
    return pga1_mul(pga1_from_re(a), b);
}

static inline pga1
pga1_lmul(pga1 a, double b)
{
    return pga1_mul(a, pga1_from_re(b));
}

static inline pga1
pga1_rdiv(double a, pga1 b)
{
    return pga1_div(pga1_from_re(a), b);
}

static inline pga1
pga1_ldiv(pga1 a, double b)
{
    return pga1_div(a, pga1_from_re(b));
}

static inline pga1
pga1_ipow(pga1 a, int64_t b)
{
    pga1 result = pga1_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return pga1_mul(a, a);
    }
    if (b == -1) {
        return pga1_inv(a);
    }
    if (b < 0) {
        a = pga1_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = pga1_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = pga1_mul(a, a);
        }
    }
    return result;
}

static inline pga1
pga1_rpow(pga1 a, pga1 b)
{
    return pga1_exp(pga1_mul(pga1_log(a), b));
}

static inline pga1
pga1_pow(pga1 a, pga1 b)
{
    return pga1_rpow(a, b);
}

static inline pga1
pga1_lpow(pga1 a, pga1 b)
{
    return pga1_exp(pga1_mul(b, pga1_log(a)));
}

static inline pga1
pga1_rpows(double a, pga1 b)
{
    return pga1_rpow(pga1_from_re(a), b);
}

static inline pga1
pga1_lpows(pga1 a, double b)
{
    return pga1_rpow(a, pga1_from_re(b));
}

static inline pga1
pga1_mix(pga1 a, pga1 b, double t)
{
    pga1 delta = pga1_mul(pga1_inv(a), b);
    return pga1_mul(a, pga1_lpows(delta, t));
}

static inline void
pga1_lmat(pga1 a, double result[ALGPY_PGA1_N*ALGPY_PGA1_N])
{
    pga1 basis = {{0}};
    pga1 product;
    int row, col;
    for (col = 0; col < ALGPY_PGA1_N; col++) {
        basis.data[col] = 1;
        product = pga1_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_PGA1_N; row++) {
            result[row * ALGPY_PGA1_N + col] = product.data[row];
        }
    }
}

static inline void
pga1_rmat(pga1 a, double result[ALGPY_PGA1_N*ALGPY_PGA1_N])
{
    pga1 basis = {{0}};
    pga1 product;
    int row, col;
    for (row = 0; row < ALGPY_PGA1_N; row++) {
        basis.data[row] = 1;
        product = pga1_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_PGA1_N; col++) {
            result[row * ALGPY_PGA1_N + col] = product.data[col];
        }
    }
}

/* Lattice API. The exterior and regressive products and the Poincare duals
 * do not depend on the metric, so they are those of exterior2, which shares
 * this coefficient order: join with identity bottom = 1, meet with identity
 * top = the pseudoscalar, rdual(e_a) join e_a = top = e_a join ldual(e_a). */
static inline pga1
pga1_bottom(void)
{
    return exterior2_bottom();
}

static inline pga1
pga1_top(void)
{
    return exterior2_top();
}

static inline pga1
pga1_join(pga1 a, pga1 b)
{
    return exterior2_join(a, b);
}

static inline pga1
pga1_meet(pga1 a, pga1 b)
{
    return exterior2_meet(a, b);
}

static inline pga1
pga1_ldual(pga1 a)
{
    return exterior2_ldual(a);
}

static inline pga1
pga1_rdual(pga1 a)
{
    return exterior2_rdual(a);
}

#endif /* ALGPY_PGA1_H */
