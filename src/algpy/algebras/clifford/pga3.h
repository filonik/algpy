/* Projective geometric algebra in 3D.
 * Generator squares: 0, 1, 1, 1.
 * 4 generators; 16 real coefficients. */
#ifndef ALGPY_PGA3_H
#define ALGPY_PGA3_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../fields/complex.h"
#include "../../fields/double.h"
#include "../../polynomials/hermite.h"
#include "pga2.h"
#include "../../vectorspaces/vec16.h"
#include "../../lattices/clifford/exterior4.h"

#define ALGPY_PGA3_N 16

typedef vec16 pga3;

static inline pga3
pga3_make(const double coeffs[ALGPY_PGA3_N])
{
    pga3 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline pga3
pga3_grade0(pga3 a)
{
    const double *x = a.data;
    return pga3_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline pga3
pga3_grade1(pga3 a)
{
    const double *x = a.data;
    return pga3_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0});
}

static inline pga3
pga3_grade2(pga3 a)
{
    const double *x = a.data;
    return pga3_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0});
}

static inline pga3
pga3_grade3(pga3 a)
{
    const double *x = a.data;
    return pga3_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0});
}

static inline pga3
pga3_grade4(pga3 a)
{
    const double *x = a.data;
    return pga3_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15]});
}

static inline pga3
pga3_from_re(double a)
{
    return pga3_make((const double[ALGPY_PGA3_N]){a});
}

static inline double
pga3_to_re(pga3 a)
{
    return a.data[0];
}

static inline pga3
pga3_zero(void)
{
    return vec16_zero();
}

static inline pga3
pga3_one(void)
{
    return pga3_from_re(1);
}

static inline pga3
pga3_pos(pga3 a)
{
    return vec16_pos(a);
}

static inline pga3
pga3_neg(pga3 a)
{
    return vec16_neg(a);
}

static inline pga3
pga3_add(pga3 a, pga3 b)
{
    return vec16_add(a, b);
}

static inline pga3
pga3_sub(pga3 a, pga3 b)
{
    return vec16_sub(a, b);
}

static inline int
pga3_equal(pga3 a, pga3 b)
{
    return vec16_all_eq(a, b);
}

static inline int
pga3_not_equal(pga3 a, pga3 b)
{
    return !pga3_equal(a, b);
}

static inline int
pga3_is_real(pga3 a)
{
    return vec16_is_real(a);
}

static inline pga3
pga3_re(pga3 a)
{
    return pga3_from_re(a.data[0]);
}

static inline pga3
pga3_im(pga3 a)
{
    return pga3_sub(a, pga3_re(a));
}

static inline pga3
pga3_mul(pga3 a, pga3 b)
{
    const double *x = a.data, *y = b.data;
    return pga3_make((const double[]){
        x[0]*y[0] + x[2]*y[2] + x[4]*y[4] - x[6]*y[6] + x[8]*y[8] - x[10]*y[10] - x[12]*y[12] - x[14]*y[14],
        x[0]*y[1] + x[1]*y[0] - x[2]*y[3] + x[3]*y[2] - x[4]*y[5] + x[5]*y[4] - x[6]*y[7] - x[7]*y[6] - x[8]*y[9] + x[9]*y[8] - x[10]*y[11] - x[11]*y[10] - x[12]*y[13] - x[13]*y[12] + x[14]*y[15] - x[15]*y[14],
        x[0]*y[2] + x[2]*y[0] - x[4]*y[6] + x[6]*y[4] - x[8]*y[10] + x[10]*y[8] - x[12]*y[14] - x[14]*y[12],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0] + x[4]*y[7] - x[5]*y[6] + x[6]*y[5] + x[7]*y[4] + x[8]*y[11] - x[9]*y[10] + x[10]*y[9] + x[11]*y[8] - x[12]*y[15] - x[13]*y[14] + x[14]*y[13] - x[15]*y[12],
        x[0]*y[4] + x[2]*y[6] + x[4]*y[0] - x[6]*y[2] - x[8]*y[12] + x[10]*y[14] + x[12]*y[8] + x[14]*y[10],
        x[0]*y[5] + x[1]*y[4] - x[2]*y[7] + x[3]*y[6] - x[4]*y[1] + x[5]*y[0] - x[6]*y[3] - x[7]*y[2] + x[8]*y[13] - x[9]*y[12] + x[10]*y[15] + x[11]*y[14] + x[12]*y[9] + x[13]*y[8] - x[14]*y[11] + x[15]*y[10],
        x[0]*y[6] + x[2]*y[4] - x[4]*y[2] + x[6]*y[0] + x[8]*y[14] - x[10]*y[12] + x[12]*y[10] + x[14]*y[8],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0] - x[8]*y[15] + x[9]*y[14] - x[10]*y[13] - x[11]*y[12] + x[12]*y[11] + x[13]*y[10] - x[14]*y[9] + x[15]*y[8],
        x[0]*y[8] + x[2]*y[10] + x[4]*y[12] - x[6]*y[14] + x[8]*y[0] - x[10]*y[2] - x[12]*y[4] - x[14]*y[6],
        x[0]*y[9] + x[1]*y[8] - x[2]*y[11] + x[3]*y[10] - x[4]*y[13] + x[5]*y[12] - x[6]*y[15] - x[7]*y[14] - x[8]*y[1] + x[9]*y[0] - x[10]*y[3] - x[11]*y[2] - x[12]*y[5] - x[13]*y[4] + x[14]*y[7] - x[15]*y[6],
        x[0]*y[10] + x[2]*y[8] - x[4]*y[14] + x[6]*y[12] - x[8]*y[2] + x[10]*y[0] - x[12]*y[6] - x[14]*y[4],
        x[0]*y[11] + x[1]*y[10] - x[2]*y[9] + x[3]*y[8] + x[4]*y[15] - x[5]*y[14] + x[6]*y[13] + x[7]*y[12] + x[8]*y[3] - x[9]*y[2] + x[10]*y[1] + x[11]*y[0] - x[12]*y[7] - x[13]*y[6] + x[14]*y[5] - x[15]*y[4],
        x[0]*y[12] + x[2]*y[14] + x[4]*y[8] - x[6]*y[10] - x[8]*y[4] + x[10]*y[6] + x[12]*y[0] + x[14]*y[2],
        x[0]*y[13] + x[1]*y[12] - x[2]*y[15] + x[3]*y[14] - x[4]*y[9] + x[5]*y[8] - x[6]*y[11] - x[7]*y[10] + x[8]*y[5] - x[9]*y[4] + x[10]*y[7] + x[11]*y[6] + x[12]*y[1] + x[13]*y[0] - x[14]*y[3] + x[15]*y[2],
        x[0]*y[14] + x[2]*y[12] - x[4]*y[10] + x[6]*y[8] + x[8]*y[6] - x[10]*y[4] + x[12]*y[2] + x[14]*y[0],
        x[0]*y[15] + x[1]*y[14] - x[2]*y[13] + x[3]*y[12] + x[4]*y[11] - x[5]*y[10] + x[6]*y[9] + x[7]*y[8] - x[8]*y[7] + x[9]*y[6] - x[10]*y[5] - x[11]*y[4] + x[12]*y[3] + x[13]*y[2] - x[14]*y[1] + x[15]*y[0]});
}

static inline pga3
pga3_square(pga3 a)
{
    const double *x = a.data;
    return pga3_make((const double[]){
        x[0]*x[0] + x[2]*x[2] + x[4]*x[4] - x[6]*x[6] + x[8]*x[8] - x[10]*x[10] - x[12]*x[12] - x[14]*x[14],
        2*x[0]*x[1] - 2*x[6]*x[7] - 2*x[10]*x[11] - 2*x[12]*x[13],
        2*x[0]*x[2] - 2*x[12]*x[14],
        2*x[0]*x[3] + 2*x[4]*x[7] + 2*x[8]*x[11] - 2*x[12]*x[15],
        2*x[0]*x[4] + 2*x[10]*x[14],
        2*x[0]*x[5] - 2*x[2]*x[7] + 2*x[8]*x[13] + 2*x[10]*x[15],
        2*x[0]*x[6] + 2*x[8]*x[14],
        2*x[0]*x[7] + 2*x[1]*x[6] - 2*x[2]*x[5] + 2*x[3]*x[4],
        2*x[0]*x[8] - 2*x[6]*x[14],
        2*x[0]*x[9] - 2*x[2]*x[11] - 2*x[4]*x[13] - 2*x[6]*x[15],
        2*x[0]*x[10] - 2*x[4]*x[14],
        2*x[0]*x[11] + 2*x[1]*x[10] - 2*x[2]*x[9] + 2*x[3]*x[8],
        2*x[0]*x[12] + 2*x[2]*x[14],
        2*x[0]*x[13] + 2*x[1]*x[12] - 2*x[4]*x[9] + 2*x[5]*x[8],
        2*x[0]*x[14] + 2*x[2]*x[12] - 2*x[4]*x[10] + 2*x[6]*x[8],
        2*x[0]*x[15] + 2*x[3]*x[12] - 2*x[5]*x[10] + 2*x[6]*x[9]});
}

static inline pga3
pga3_conjugate(pga3 a)
{
    const double *x = a.data;
    return pga3_make((const double[]){x[0], -x[1], -x[2], -x[3], -x[4], -x[5], -x[6], x[7], -x[8], -x[9], -x[10], x[11], -x[12], x[13], x[14], x[15]});
}

static inline void
pga3_basis(pga3 result[16])
{
    int i, j;
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
pga3_gens(pga3 result[4])
{
    result[0] = pga3_make((const double[ALGPY_PGA3_N]){[1] = 1});
    result[1] = pga3_make((const double[ALGPY_PGA3_N]){[2] = 1});
    result[2] = pga3_make((const double[ALGPY_PGA3_N]){[4] = 1});
    result[3] = pga3_make((const double[ALGPY_PGA3_N]){[8] = 1});
}

static inline void
pga3_coeffs(pga3 result[256])
{
    pga3 basis[16];
    int i, j;
    pga3_basis(basis);
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            result[i*16+j] = pga3_mul(basis[i], basis[j]);
        }
    }
}

static inline double
pga3_scalar_product(pga3 a, pga3 b)
{
    const double *x = a.data, *y = b.data;
    return x[0]*y[0] + x[2]*y[2] + x[4]*y[4] - x[6]*y[6] + x[8]*y[8] - x[10]*y[10] - x[12]*y[12] - x[14]*y[14];
}

/* Centred scalar traces determine the generic quartic, including
 * multiplicities for elements with a lower-degree element relation. */
static inline void
pga3_minpoly(pga3 a, double result[5])
{
    double w = a.data[0], s2, s3, s4;
    pga3 b = a, b2;
    b.data[0] = 0;
    b2 = pga3_square(b);
    s2 = b2.data[0];
    s3 = pga3_scalar_product(b2, b);
    s4 = pga3_scalar_product(b2, b2);
    result[0] = w*w*w*w - 2*w*w*s2 + 4*w*s3/3 + 2*s2*s2 - s4;
    result[1] = -4*w*w*w + 4*w*s2 - 4*s3/3;
    result[2] = 6*w*w - 2*s2;
    result[3] = -4*w;
    result[4] = 1;
}

/* The null radical drops out. The Euclidean Cl(3,0) quotient has a
 * complex determinant d0 + d7*I, so the quartic norm is d0*d0 + d7*d7. */
static inline double
pga3_norm4(pga3 a)
{
    const double *x = a.data;
    double d0 = x[0]*x[0] - x[14]*x[14] - x[2]*x[2] - x[4]*x[4]
        - x[8]*x[8] + x[6]*x[6] + x[10]*x[10] + x[12]*x[12];
    double d7 = 2*(x[0]*x[14] - x[2]*x[12] + x[4]*x[10] - x[8]*x[6]);
    return d0*d0 + d7*d7;
}

static inline double
pga3_norm(pga3 a)
{
    const double *x = a.data;
    double d0 = x[0]*x[0] - x[14]*x[14] - x[2]*x[2] - x[4]*x[4]
        - x[8]*x[8] + x[6]*x[6] + x[10]*x[10] + x[12]*x[12];
    double d7 = 2*(x[0]*x[14] - x[2]*x[12] + x[4]*x[10] - x[8]*x[6]);
    return sqrt(hypot(d0, d7));
}

static inline pga3
pga3_inv(pga3 a)
{
    double p[5];
    pga3 a2, a3, r;
    pga3_minpoly(a, p);
    a2 = pga3_square(a);
    a3 = pga3_mul(a2, a);
    r = vec16_add(vec16_add(a3, vec16_scale(a2, p[3])),
                  vec16_add(vec16_scale(a, p[2]), pga3_from_re(p[1])));
    return vec16_div_scalar(vec16_neg(r), p[0]);
}

static inline void
pga3_eivals(pga3 a, complex result[4])
{
    double p[5];
    pga3_minpoly(a, p);
    poly4_coeffs_to_roots(p, result);
}

static inline pga3
pga3_div(pga3 a, pga3 b)
{
    return pga3_mul(a, pga3_inv(b));
}

static inline pga3
pga3_real_slice(complex value)
{
    pga3 result = pga3_from_re(value.data[0]);
    result.data[6] = value.data[1];
    return result;
}

static inline int
pga3_unused_generator(pga3 a)
{
    int generator, i;
    for (generator = 1; generator < 4; generator++) {
        int used = 0;
        for (i = 1; i < 16; i++) {
            if ((i & (1 << generator)) && a.data[i] != 0) {
                used = 1;
                break;
            }
        }
        if (!used) {
            return generator;
        }
    }
    return -1;
}

/* All regular generators have square +1. An unused one can be removed
 * without changing the metric or reordering the special null generator. */
static inline pga3
pga3_apply_subalgebra(pga3 a, int generator,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    pga2 suba, value;
    pga3 result = {{0}};
    int low = (1 << generator) - 1, i;
    for (i = 0; i < 8; i++) {
        int mask = (i & low) | ((i & ~low) << 1);
        suba.data[i] = a.data[mask];
    }
    value = pga2_apply(suba, f, deriv1, deriv2, deriv3, cf, cderiv1);
    if (isnan(value.data[0])) {
        return pga3_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
    }
    for (i = 0; i < 8; i++) {
        int mask = (i & low) | ((i & ~low) << 1);
        result.data[mask] = value.data[i];
    }
    return result;
}

static inline int
pga3_quadratic_relation(pga3 b, pga3 b2, double *linear)
{
    int i, pivot = 1;
    double bound, size = 0;
    for (i = 2; i < 16; i++) {
        if (fabs(b.data[i]) > fabs(b.data[pivot])) {
            pivot = i;
        }
    }
    *linear = b2.data[pivot]/b.data[pivot];
    if (!isfinite(*linear)) {
        return 0;
    }
    for (i = 0; i < 16; i++) {
        size = fmax(size, fabs(b2.data[i]));
    }
    bound = 32*DBL_EPSILON*(size + fabs(*linear*b.data[pivot]));
    for (i = 1; i < 16; i++) {
        double value = *linear*b.data[i];
        if (fabs(b2.data[i] - value) > bound) {
            return 0;
        }
    }
    return 1;
}

static inline pga3
pga3_apply_quadratic(double w, pga3 b, double constant, double linear,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    double shift = linear/2;
    if (hermite_quadratic_apply(w + shift, constant + shift*shift, 0,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        pga3 result = vec16_scale(b, value.odd);
        result.data[0] = value.even - shift*value.odd;
        return result;
    }
    return pga3_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                          NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
}

static inline int
pga3_cubic_relation(pga3 b, pga3 b2, pga3 b3, double p[4])
{
    int i, first = 1, second = 1;
    double ratio, largest = 0, linear, quadratic, bound = 0;
    for (i = 2; i < 16; i++) {
        if (fabs(b.data[i]) > fabs(b.data[first])) {
            first = i;
        }
    }
    ratio = b2.data[first]/b.data[first];
    for (i = 1; i < 16; i++) {
        double residual = fabs(b2.data[i] - ratio*b.data[i]);
        if (residual > largest) {
            largest = residual;
            second = i;
        }
    }
    if (largest == 0) {
        return 0;
    }
    quadratic = (b3.data[second] - b3.data[first]/b.data[first]*b.data[second])/
        (b2.data[second] - ratio*b.data[second]);
    linear = (b3.data[first] - quadratic*b2.data[first])/b.data[first];
    if (!isfinite(linear) || !isfinite(quadratic)) {
        return 0;
    }
    for (i = 0; i < 16; i++) {
        bound = fmax(bound, fabs(b3.data[i]) + fabs(quadratic*b2.data[i]) + fabs(linear*b.data[i]));
    }
    bound *= 64*DBL_EPSILON;
    for (i = 1; i < 16; i++) {
        double q = quadratic*b2.data[i], l = linear*b.data[i];
        if (fabs(b3.data[i] - q - l) > bound) {
            return 0;
        }
    }
    p[0] = quadratic*b2.data[0] - b3.data[0];
    p[1] = -linear;
    p[2] = -quadratic;
    p[3] = 1;
    return 1;
}

static inline pga3
pga3_apply(pga3 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double w = a.data[0], linear, coefficients[4], cubic[4];
    pga3 b = a, b2, b3, result;
    int generator, ok, i;
    if (vec16_is_real(a)) {
        return pga3_real_slice(cf(complex_make((const double[]){w, 0})));
    }
    generator = pga3_unused_generator(a);
    if (generator >= 0) {
        return pga3_apply_subalgebra(a, generator,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    b.data[0] = 0;
    b2 = pga3_square(b);
    if (pga3_quadratic_relation(b, b2, &linear)) {
        return pga3_apply_quadratic(w, b, b2.data[0], linear,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    b3 = pga3_mul(b2, b);
    if (pga3_cubic_relation(b, b2, b3, cubic)) {
        complex roots[4];
        poly3_coeffs_to_roots(cubic, roots);
        ok = hermite_interpolate(3, w, roots, f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients);
    } else {
        double s2 = b2.data[0], s3 = b3.data[0];
        double s4 = pga3_scalar_product(b2, b2);
        double quartic[5] = {2*s2*s2 - s4, -4*s3/3, -2*s2, 0, 1};
        ok = poly4_apply(quartic, w, f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients);
    }
    if (!ok) {
        return pga3_make((const double[]){NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                              NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
    }
    for (i = 0; i < 16; i++) {
        result.data[i] = coefficients[1]*b.data[i] + coefficients[2]*b2.data[i] + coefficients[3]*b3.data[i];
    }
    result.data[0] += coefficients[0];
    return result;
}

static inline pga3
pga3_exp(pga3 a)
{
    return pga3_apply(a, double_exp, double_exp_deriv1,
        double_exp_deriv2, double_exp_deriv3, complex_exp, complex_exp_deriv1);
}

static inline pga3
pga3_log(pga3 a)
{
    return pga3_apply(a, double_log, double_log_deriv1,
        double_log_deriv2, double_log_deriv3, complex_log, complex_log_deriv1);
}

static inline pga3
pga3_sqrt(pga3 a)
{
    return pga3_apply(a, double_sqrt, double_sqrt_deriv1,
        double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt, complex_sqrt_deriv1);
}

static inline pga3
pga3_cos(pga3 a)
{
    return pga3_apply(a, double_cos, double_cos_deriv1,
        double_cos_deriv2, double_cos_deriv3, complex_cos, complex_cos_deriv1);
}

static inline pga3
pga3_sin(pga3 a)
{
    return pga3_apply(a, double_sin, double_sin_deriv1,
        double_sin_deriv2, double_sin_deriv3, complex_sin, complex_sin_deriv1);
}

static inline pga3
pga3_tan(pga3 a)
{
    return pga3_apply(a, double_tan, double_tan_deriv1,
        double_tan_deriv2, double_tan_deriv3, complex_tan, complex_tan_deriv1);
}

static inline pga3
pga3_acos(pga3 a)
{
    return pga3_apply(a, double_acos, double_acos_deriv1,
        double_acos_deriv2, double_acos_deriv3, complex_acos, complex_acos_deriv1);
}

static inline pga3
pga3_asin(pga3 a)
{
    return pga3_apply(a, double_asin, double_asin_deriv1,
        double_asin_deriv2, double_asin_deriv3, complex_asin, complex_asin_deriv1);
}

static inline pga3
pga3_atan(pga3 a)
{
    return pga3_apply(a, double_atan, double_atan_deriv1,
        double_atan_deriv2, double_atan_deriv3, complex_atan, complex_atan_deriv1);
}

static inline pga3
pga3_cosh(pga3 a)
{
    return pga3_apply(a, double_cosh, double_cosh_deriv1,
        double_cosh_deriv2, double_cosh_deriv3, complex_cosh, complex_cosh_deriv1);
}

static inline pga3
pga3_sinh(pga3 a)
{
    return pga3_apply(a, double_sinh, double_sinh_deriv1,
        double_sinh_deriv2, double_sinh_deriv3, complex_sinh, complex_sinh_deriv1);
}

static inline pga3
pga3_tanh(pga3 a)
{
    return pga3_apply(a, double_tanh, double_tanh_deriv1,
        double_tanh_deriv2, double_tanh_deriv3, complex_tanh, complex_tanh_deriv1);
}

static inline pga3
pga3_acosh(pga3 a)
{
    return pga3_apply(a, double_acosh, double_acosh_deriv1,
        double_acosh_deriv2, double_acosh_deriv3, complex_acosh, complex_acosh_deriv1);
}

static inline pga3
pga3_asinh(pga3 a)
{
    return pga3_apply(a, double_asinh, double_asinh_deriv1,
        double_asinh_deriv2, double_asinh_deriv3, complex_asinh, complex_asinh_deriv1);
}

static inline pga3
pga3_atanh(pga3 a)
{
    return pga3_apply(a, double_atanh, double_atanh_deriv1,
        double_atanh_deriv2, double_atanh_deriv3, complex_atanh, complex_atanh_deriv1);
}

static inline pga3
pga3_expm1(pga3 a)
{
    return pga3_apply(a, double_expm1, double_expm1_deriv1,
        double_expm1_deriv2, double_expm1_deriv3, complex_expm1, complex_expm1_deriv1);
}

static inline pga3
pga3_log1p(pga3 a)
{
    return pga3_apply(a, double_log1p, double_log1p_deriv1,
        double_log1p_deriv2, double_log1p_deriv3, complex_log1p, complex_log1p_deriv1);
}

static inline pga3
pga3_radd(double a, pga3 b)
{
    return pga3_add(pga3_from_re(a), b);
}

static inline pga3
pga3_ladd(pga3 a, double b)
{
    return pga3_add(a, pga3_from_re(b));
}

static inline pga3
pga3_rsub(double a, pga3 b)
{
    return pga3_sub(pga3_from_re(a), b);
}

static inline pga3
pga3_lsub(pga3 a, double b)
{
    return pga3_sub(a, pga3_from_re(b));
}

static inline pga3
pga3_rmul(double a, pga3 b)
{
    return pga3_mul(pga3_from_re(a), b);
}

static inline pga3
pga3_lmul(pga3 a, double b)
{
    return pga3_mul(a, pga3_from_re(b));
}

static inline pga3
pga3_rdiv(double a, pga3 b)
{
    return pga3_div(pga3_from_re(a), b);
}

static inline pga3
pga3_ldiv(pga3 a, double b)
{
    return pga3_div(a, pga3_from_re(b));
}

static inline pga3
pga3_ipow(pga3 a, int64_t b)
{
    pga3 result = pga3_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return pga3_mul(a, a);
    }
    if (b == -1) {
        return pga3_inv(a);
    }
    if (b < 0) {
        a = pga3_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = pga3_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = pga3_mul(a, a);
        }
    }
    return result;
}

static inline pga3
pga3_rpow(pga3 a, pga3 b)
{
    return pga3_exp(pga3_mul(pga3_log(a), b));
}

static inline pga3
pga3_pow(pga3 a, pga3 b)
{
    return pga3_rpow(a, b);
}

static inline pga3
pga3_lpow(pga3 a, pga3 b)
{
    return pga3_exp(pga3_mul(b, pga3_log(a)));
}

static inline pga3
pga3_rpows(double a, pga3 b)
{
    return pga3_rpow(pga3_from_re(a), b);
}

static inline pga3
pga3_lpows(pga3 a, double b)
{
    return pga3_rpow(a, pga3_from_re(b));
}

static inline pga3
pga3_mix(pga3 a, pga3 b, double t)
{
    pga3 delta = pga3_mul(pga3_inv(a), b);
    return pga3_mul(a, pga3_lpows(delta, t));
}

static inline void
pga3_lmat(pga3 a, double result[ALGPY_PGA3_N*ALGPY_PGA3_N])
{
    pga3 basis = {{0}};
    pga3 product;
    int row, col;
    for (col = 0; col < ALGPY_PGA3_N; col++) {
        basis.data[col] = 1;
        product = pga3_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_PGA3_N; row++) {
            result[row * ALGPY_PGA3_N + col] = product.data[row];
        }
    }
}

static inline void
pga3_rmat(pga3 a, double result[ALGPY_PGA3_N*ALGPY_PGA3_N])
{
    pga3 basis = {{0}};
    pga3 product;
    int row, col;
    for (row = 0; row < ALGPY_PGA3_N; row++) {
        basis.data[row] = 1;
        product = pga3_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_PGA3_N; col++) {
            result[row * ALGPY_PGA3_N + col] = product.data[col];
        }
    }
}

/* Lattice API. The exterior and regressive products and the Poincare duals
 * do not depend on the metric, so they are those of exterior4, which shares
 * this coefficient order: join with identity bottom = 1, meet with identity
 * top = the pseudoscalar, rdual(e_a) join e_a = top = e_a join ldual(e_a). */
static inline pga3
pga3_bottom(void)
{
    return exterior4_bottom();
}

static inline pga3
pga3_top(void)
{
    return exterior4_top();
}

static inline pga3
pga3_join(pga3 a, pga3 b)
{
    return exterior4_join(a, b);
}

static inline pga3
pga3_meet(pga3 a, pga3 b)
{
    return exterior4_meet(a, b);
}

static inline pga3
pga3_ldual(pga3 a)
{
    return exterior4_ldual(a);
}

static inline pga3
pga3_rdual(pga3 a)
{
    return exterior4_rdual(a);
}

#endif /* ALGPY_PGA3_H */
