#ifndef ALGPY_HERMITE_QUADRATIC_H
#define ALGPY_HERMITE_QUADRATIC_H

#include "../fields/double.h"
#include "../fields/complex.h"

/* Distance to a native function's nearest cut or pole. A negative result
 * denotes a custom callback, for which no length scale is assumed. */
static inline double
hermite_function_distance(complex z, double (*f)(double))
{
    double x = z.data[0], y = z.data[1];
    if (f == double_log || f == double_sqrt || f == double_inv) {
        return x <= 0 && f != double_inv ? fabs(y) : hypot(x, y);
    }
    if (f == double_log1p) {
        return x <= -1 ? fabs(y) : hypot(x + 1, y);
    }
    if (f == double_asin || f == double_acos || f == double_atanh) {
        return fabs(x) >= 1 ? fabs(y) : fmin(hypot(x - 1, y), hypot(x + 1, y));
    }
    if (f == double_acosh) {
        return x <= 1 ? fabs(y) : hypot(x - 1, y);
    }
    if (f == double_asinh || f == double_atan) {
        return fabs(y) >= 1 ? fabs(x) : fmin(hypot(x, y - 1), hypot(x, y + 1));
    }
    if (f == double_tan) {
        return hypot(remainder(x - M_PI/2, M_PI), y);
    }
    if (f == double_tanh) {
        return hypot(x, remainder(y - M_PI/2, M_PI));
    }
    if (f == double_exp || f == double_expm1 || f == double_sin ||
        f == double_cos || f == double_sinh || f == double_cosh || f == double_square) {
        return 1;
    }
    return -1;
}

/* f(w+u) = even + odd*u for u*u=q. Subscripts dw and dq denote partial
 * derivatives, needed only by dual centres and their interpolation formulas. */
typedef struct {
    double even, odd;
    double even_dw, odd_dw, odd_dq;
} hermite_quadratic;

static inline int
hermite_quadratic_is_finite(const hermite_quadratic *value, int with_derivatives)
{
    return isfinite(value->even) && isfinite(value->odd) &&
        (!with_derivatives || (isfinite(value->even_dw) &&
            isfinite(value->odd_dw) && isfinite(value->odd_dq)));
}

/* Coalescing real or conjugate roots. Derivative samples avoid
 * subtracting nearly equal values and retain the terms through q^3.
 * Return 0 if the required real derivatives are not finite. */
static inline int
hermite_quadratic_near(double w, double q, double radius, int with_derivatives,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double), hermite_quadratic *result)
{
    double d1, d2, d3;
    d1 = deriv1(w);
    d2 = deriv2(w);
    d3 = deriv3(w);
    if (!(isfinite(f(w)) && isfinite(d1) && isfinite(d2) && isfinite(d3))) {
        return 0;
    }
    /* Integrals of f', f'' and f''' determine B, B_w and B_q.
     * These bounded weights reproduce their first four even Taylor
     * terms using the centre and paired samples at r, 2r and 3r. */
    static const double even_plus[4] = {586./945, 167./420, -2./105, 1./756};
    static const double even_minus[4] = {1451./945, -43./70, 3./35, -13./1890};
    static const double mean_plus[4] = {12067./30240, 2171./20160, -73./10080, 31./60480};
    static const double mean_minus[4] = {6319./10080, -2899./20160, 37./2016, -29./20160};
    static const double slope_plus[4] = {359./2835, 23./540, -13./4725, 11./56700};
    static const double slope_minus[4] = {247./1134, -439./7560, 71./9450, -67./113400};
    const double *even = q > 0 ? even_plus : even_minus;
    const double *mean = q > 0 ? mean_plus : mean_minus;
    const double *slope = q > 0 ? slope_plus : slope_minus;
    double average = mean[0]*d2;
    int i;
    result->odd = even[0]*d1;
    result->even_dw = mean[0]*d3;
    result->odd_dw = even[0]*d2;
    result->odd_dq = slope[0]*d3;
    for (i = 1; i <= 3; i++) {
        double offset = i*radius;
        double first = deriv1(w + offset)/2 + deriv1(w - offset)/2;
        double second = deriv2(w + offset)/2 + deriv2(w - offset)/2;
        result->odd += even[i]*first;
        average += mean[i]*second;
        if (with_derivatives) {
            double third = deriv3(w + offset)/2 + deriv3(w - offset)/2;
            result->even_dw += mean[i]*third;
            result->odd_dw += even[i]*second;
            result->odd_dq += slope[i]*third;
        }
    }
    result->even = f(w) + q*average;
    result->even_dw = d1 + q*result->even_dw;
    return 1;
}

/* Real quadratic values, with optional w/q derivatives for dual centres.
 * even_dq = odd_dw/2, so only odd_dq needs a separate field. */
static inline int
hermite_quadratic_apply(double w, double q, int with_derivatives,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), hermite_quadratic *result)
{
    double radius = sqrt(fabs(q));
    double distance = hermite_function_distance(complex_make(w, 0), f);
    if (q != 0 && radius < 0.015*distance &&
        hermite_quadratic_near(w, q, radius, with_derivatives,
            f, deriv1, deriv2, deriv3, result)) {
        return hermite_quadratic_is_finite(result, with_derivatives);
    }
    if (q == 0) {
        result->even = f(w);
        result->odd = deriv1(w);
        result->even_dw = result->odd;
        result->odd_dw = with_derivatives ? deriv2(w) : 0;
        result->odd_dq = with_derivatives ? deriv3(w)/6 : 0;
    } else if (q > 0) {
        double p = f(w + radius), m = f(w - radius);
        result->even = p/2 + m/2;
        result->odd = (p/2 - m/2)/radius;
        if (with_derivatives) {
            double dp = deriv1(w + radius), dm = deriv1(w - radius);
            result->even_dw = dp/2 + dm/2;
            result->odd_dw = (dp/2 - dm/2)/radius;
            result->odd_dq = (result->even_dw - result->odd)/(2*q);
        }
    } else if (q < 0) {
        complex value = cf(complex_make(w, radius));
        result->even = value.data[0];
        result->odd = value.data[1]/radius;
        if (with_derivatives) {
            complex derivative = cderiv1(complex_make(w, radius));
            result->even_dw = derivative.data[0];
            result->odd_dw = derivative.data[1]/radius;
            result->odd_dq = (result->even_dw - result->odd)/(2*q);
        }
    } else {
        return 0;
    }
    return hermite_quadratic_is_finite(result, with_derivatives);
}

/* Integral of f' along z+[-root,root], divided by the segment length.
 * Four-point Gauss-Legendre quadrature avoids subtracting close values. */
static inline complex
hermite_first_difference(complex z, complex root, complex (*deriv1)(complex))
{
    complex offset = vec2_scale(root, 0.3399810435848562648);
    complex p = deriv1(vec2_add(z, offset));
    complex m = deriv1(vec2_sub(z, offset));
    complex result = vec2_scale(vec2_add(p, m), 0.3260725774312730713);
    offset = vec2_scale(root, 0.8611363115940525752);
    p = deriv1(vec2_add(z, offset));
    m = deriv1(vec2_sub(z, offset));
    return vec2_add(result, vec2_scale(vec2_add(p, m), 0.1739274225687269287));
}

/* Coalescing complex roots: integrate f' instead of forming a divided
 * difference. The even formula retains terms through q^3; four-point
 * Gauss-Legendre quadrature gives the same accuracy for the odd coefficient. */
static inline void
hermite_complex_quadratic_near(complex z, complex root,
    complex (*cf)(complex), complex (*cderiv1)(complex), complex *even, complex *odd)
{
    complex offset = vec2_scale(root, 1/sqrt(3));
    complex p = cderiv1(vec2_add(z, offset));
    complex m = cderiv1(vec2_sub(z, offset));
    complex difference = vec2_scale(vec2_sub(p, m), 3*sqrt(3)/16);
    p = cderiv1(vec2_add(z, root));
    m = cderiv1(vec2_sub(z, root));
    difference = vec2_add(difference, vec2_scale(vec2_sub(p, m), 1./16));
    *even = vec2_add(cf(z), complex_mul(root, difference));
    *odd = hermite_first_difference(z, root, cderiv1);
}

#endif /* ALGPY_HERMITE_QUADRATIC_H */
