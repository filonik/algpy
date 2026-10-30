#ifndef ALGPY_CLIFFORD3_H
#define ALGPY_CLIFFORD3_H

#define ALGPY_CLIFFORD3_N 8

#include <stdint.h>
#include <math.h>

#include "../../../fields/complex.h"
#include "../../../fields/double.h"
#include "../../../polynomials/hermite_quadratic.h"
#include "../../../vectorspaces/vec3.h"
#include "../../../vectorspaces/vec8.h"

typedef vec8 clifford3;

static inline clifford3
clifford3_make(double x, double x0, double x1, double x01,
                    double x2, double x02, double x12, double x012)
{
    clifford3 r = {{x, x0, x1, x01, x2, x02, x12, x012}};
    return r;
}

static inline clifford3
clifford3_grade0(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(x[0], 0, 0, 0, 0, 0, 0, 0);
}

static inline clifford3
clifford3_grade1(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(0, x[1], x[2], 0, x[4], 0, 0, 0);
}

static inline clifford3
clifford3_grade2(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(0, 0, 0, x[3], 0, x[5], x[6], 0);
}

static inline clifford3
clifford3_grade3(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(0, 0, 0, 0, 0, 0, 0, x[7]);
}

static inline clifford3
clifford3_from_re(double a)
{
    return clifford3_make(a, 0, 0, 0, 0, 0, 0, 0);
}

static inline double
clifford3_to_re(clifford3 a)
{
    return a.data[0];
}

static inline clifford3
clifford3_one(void)
{
    return clifford3_from_re(1);
}

static inline clifford3
clifford3_radd(double a, clifford3 b) {
  return vec8_add(clifford3_from_re(a), b);
}

static inline clifford3
clifford3_ladd(clifford3 a, double b) {
  return vec8_add(a, clifford3_from_re(b));
}

static inline clifford3
clifford3_rsub(double a, clifford3 b) {
  return vec8_sub(clifford3_from_re(a), b);
}

static inline clifford3
clifford3_lsub(clifford3 a, double b) {
  return vec8_sub(a, clifford3_from_re(b));
}

static inline clifford3
clifford3_mul(vec3 k, clifford3 a, clifford3 b)
{
    clifford3 r = {{0}};
    int i, j, p, q;
    for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD3_N; j++) {
            double factor = 1;
            int swaps = 0;
            for (p = 0; p < 3; p++) {
                if (!(j & (1 << p))) {
                    continue;
                }
                for (q = p + 1; q < 3; q++) {
                    if (i & (1 << q)) {
                        swaps++;
                    }
                }
                if (i & (1 << p)) {
                    factor *= k.data[p];
                }
            }
            if (swaps & 1) {
                factor = -factor;
            }
            r.data[i ^ j] += factor * a.data[i] * b.data[j];
        }
    }
    return r;
}

/* Clifford conjugate: grade 1 and grade 2 components change sign. */
static inline clifford3
clifford3_conjugate(clifford3 a)
{
    const double *x = a.data;
    return clifford3_make(x[0], -x[1], -x[2], -x[3], -x[4], -x[5], -x[6], x[7]);
}

/* Product with the conjugate: lies in the centre span{1, e0*e1*e2}. */
static inline clifford3
clifford3_conjugate_product(vec3 k, clifford3 a)
{
    return clifford3_mul(k, a, clifford3_conjugate(a));
}

/* Norm (degree 4), the constant coefficient of the minimal polynomial. */
static inline double
clifford3_norm4(vec3 k, clifford3 a)
{
    clifford3 d = clifford3_conjugate_product(k, a);
    double s = -k.data[0] * k.data[1] * k.data[2];
    return d.data[0] * d.data[0] - s * d.data[7] * d.data[7];
}

static inline complex
clifford3_norm(vec3 k, clifford3 a)
{
    return complex_root_real(clifford3_norm4(k, a), 4);
}

static inline clifford3
clifford3_inv(vec3 k, clifford3 a)
{
    clifford3 d = clifford3_conjugate_product(k, a);
    double s = -k.data[0] * k.data[1] * k.data[2];
    double n = d.data[0] * d.data[0] - s * d.data[7] * d.data[7];
    clifford3 dinv = clifford3_make(d.data[0], 0, 0, 0, 0, 0, 0, -d.data[7]);
    return vec8_div_scalar(
        clifford3_mul(k, clifford3_conjugate(a), dinv), n);
}

static inline clifford3
clifford3_div(vec3 k, clifford3 a, clifford3 b)
{
    return clifford3_mul(k, a, clifford3_inv(k, b));
}

static inline void
clifford3_basis(clifford3 result[ALGPY_CLIFFORD3_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
        double components[ALGPY_CLIFFORD3_N] = {0};
        components[i] = 1;
        for (j = 0; j < ALGPY_CLIFFORD3_N; j++) {
            result[i].data[j] = components[j];
        }
    }
}

static inline void
clifford3_coeffs(
    vec3 k,
    clifford3 result[ALGPY_CLIFFORD3_N * ALGPY_CLIFFORD3_N])
{
    clifford3 basis[ALGPY_CLIFFORD3_N];
    int i, j;
    clifford3_basis(basis);
    for (i = 0; i < ALGPY_CLIFFORD3_N; i++) {
        for (j = 0; j < ALGPY_CLIFFORD3_N; j++) {
            result[i * ALGPY_CLIFFORD3_N + j] =
                clifford3_mul(k, basis[i], basis[j]);
        }
    }
}

static inline void
clifford3_minpoly(vec3 k, clifford3 a, double result[5])
{
    const int n = 4;
    clifford3 uk = a;
    int i;

    result[n] = 1;
    for (i = 1; i <= n; i++) {
        double ck = ((double)n / i) * clifford3_to_re(uk);
        result[n - i] = -ck;

        if (i < n) {
            clifford3 adj = clifford3_lsub(uk, ck);
            uk = clifford3_mul(k, a, adj);
        }
    }
}

static inline void
clifford3_eivals(vec3 k, clifford3 a, complex result[4])
{
    /*
     * The centre is span{1, z} with z = e0*e1*e2, z*z = s. Each element
     * satisfies the central quadratic t^2 - tau*t + delta, where
     * tau = 2*(x + x012*z) and delta = a*conj(a).
     */
    clifford3 d = clifford3_conjugate_product(k, a);
    double s = -k.data[0] * k.data[1] * k.data[2];
    complex w = complex_sqrt_real(s);
    int m;
    for (m = 0; m < 2; m++) {
        double sign = m == 0 ? -1 : 1;
        complex tau = complex_make(
            2 * (a.data[0] + sign * a.data[7] * w.data[0]),
            2 * sign * a.data[7] * w.data[1]);
        complex delta = complex_make(
            d.data[0] + sign * d.data[7] * w.data[0],
            sign * d.data[7] * w.data[1]);
        complex tau2 = complex_mul(tau, tau);
        complex root = complex_sqrt(complex_make(
            tau2.data[0] - 4 * delta.data[0],
            tau2.data[1] - 4 * delta.data[1]));
        result[2 * m] = complex_make((tau.data[0] - root.data[0]) / 2,
                                      (tau.data[1] - root.data[1]) / 2);
        result[2 * m + 1] = complex_make((tau.data[0] + root.data[0]) / 2,
                                          (tau.data[1] + root.data[1]) / 2);
    }
}

/* Common algebra API. */
static inline clifford3
clifford3_zero(void)
{
    return vec8_zero();
}

static inline int
clifford3_is_real(clifford3 a)
{
    return vec8_is_real(a);
}

static inline clifford3
clifford3_pos(clifford3 a)
{
    return vec8_pos(a);
}

static inline clifford3
clifford3_neg(clifford3 a)
{
    return vec8_neg(a);
}

static inline clifford3
clifford3_add(clifford3 a, clifford3 b)
{
    return vec8_add(a, b);
}

static inline clifford3
clifford3_sub(clifford3 a, clifford3 b)
{
    return vec8_sub(a, b);
}

static inline int
clifford3_equal(clifford3 a, clifford3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
clifford3_not_equal(clifford3 a, clifford3 b)
{
    return !clifford3_equal(a, b);
}

static inline clifford3
clifford3_re(clifford3 a)
{
    return clifford3_from_re(a.data[0]);
}

static inline clifford3
clifford3_im(clifford3 a)
{
    return clifford3_sub(a, clifford3_re(a));
}

/* A + B*u, where u contains the vector and bivector parts of a, and
 * A = a0 + a7*I, B = b0 + b7*I, I = e0*e1*e2. */
static inline clifford3
clifford3_center_slice(vec3 k, clifford3 a,
    double a0, double a7, double b0, double b7)
{
    const double *x = a.data;
    double k0 = k.data[0], k1 = k.data[1], k2 = k.data[2];
    return clifford3_make(
        a0,
        b0*x[1] - b7*k1*k2*x[6],
        b0*x[2] + b7*k0*k2*x[5],
        b0*x[3] + b7*k2*x[4],
        b0*x[4] - b7*k0*k1*x[3],
        b0*x[5] - b7*k1*x[2],
        b0*x[6] + b7*k0*x[1],
        a7);
}

/* Real scalar inputs use the first available negative-square basis blade
 * for the same principal complex continuation as the Cayley-Dickson family. */
static inline clifford3
clifford3_real_slice(vec3 k, complex value)
{
    clifford3 result = clifford3_from_re(value.data[0]);
    double squares[7] = {k.data[0], k.data[1], -k.data[0]*k.data[1],
        k.data[2], -k.data[0]*k.data[2], -k.data[1]*k.data[2],
        -k.data[0]*k.data[1]*k.data[2]};
    int i;
    if (value.data[1] == 0) {
        return result;
    }
    for (i = 0; i < 7; i++) {
        if (squares[i] < 0) {
            result.data[i + 1] = value.data[1]/sqrt(-squares[i]);
            return result;
        }
    }
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* f(a) = (even_re + even_i*I) + (odd_re + odd_i*I)*u. */
typedef struct {
    double even_re, even_i, odd_re, odd_i;
} clifford3_center_coefficients;

static inline int
clifford3_is_central(clifford3 a)
{
    const double *x = a.data;
    return x[1] == 0 && x[2] == 0 && x[3] == 0 &&
        x[4] == 0 && x[5] == 0 && x[6] == 0;
}

/* Detect coalescing centre channels using the nearest cut or pole. */
static inline int
clifford3_center_is_near_dual(double w, double t, double q0, double q7, double s,
    double (*f)(double))
{
    double radius, distance;
    if (s == 0 || !isfinite(s)) {
        return 0;
    }
    radius = sqrt(fabs(q0));
    if (q0 >= 0) {
        distance = fmin(hermite_function_distance(complex_make(w + radius, 0), f),
            hermite_function_distance(complex_make(w - radius, 0), f));
    } else {
        distance = hermite_function_distance(complex_make(w, radius), f);
    }
    return distance > 0 && sqrt(fabs(s))*(fabs(t)/distance +
        fabs(q7)/distance/distance) < 0.001 && isfinite(f(w));
}

/* Near I^2=0, interpolate centre coefficients through s^3 from real-axis
 * derivatives. Return 0 if an interpolation point is outside the domain. */
static inline int
clifford3_near_dual_coefficients(double w, double t, double q0, double q7, double s,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), clifford3_center_coefficients *result)
{
    double radius = sqrt(fabs(s));
    hermite_quadratic base = {0};
    static const double odd_plus[4] = {586./945, 167./420, -2./105, 1./756};
    static const double odd_minus[4] = {1451./945, -43./70, 3./35, -13./1890};
    static const double even_plus[3] = {179./576, -13./360, 11./2880};
    static const double even_minus[3] = {-257./576, 43./360, -41./2880};
    const double *odd = s > 0 ? odd_plus : odd_minus;
    const double *even = s > 0 ? even_plus : even_minus;
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
        double offset = i*radius;
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
        result->even_re += radius*even[i - 1]*(ap - am);
        result->even_i += odd[i]*(ap/2 + am/2);
        result->odd_re += radius*even[i - 1]*(bp - bm);
        result->odd_i += odd[i]*(bp/2 + bm/2);
    }
    return 1;
}

/* Elements in span{1,I} need no vector/bivector reconstruction. */
static inline clifford3
clifford3_apply_central(clifford3 a, double s,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double w = a.data[0], t = a.data[7];
    if (s < 0) {
        double scale = sqrt(-s);
        complex value = cf(complex_make(w, scale*t));
        return clifford3_make(value.data[0], 0, 0, 0, 0, 0, 0, value.data[1]/scale);
    }
    if (s > 0) {
        hermite_quadratic value = {0};
        if (hermite_quadratic_apply(w, s*t*t, 0,
                f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
            return clifford3_make(value.even, 0, 0, 0, 0, 0, 0, t*value.odd);
        }
    } else if (s == 0) {
        double value = f(w), derivative = deriv1(w);
        if (isfinite(value) && isfinite(derivative)) {
            return clifford3_make(value, 0, 0, 0, 0, 0, 0, t*derivative);
        }
    }
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* Scalar centre and scalar u^2: the same real quadratic calculus as
 * Cayley-Dickson, including square-zero nilpotents. */
static inline clifford3
clifford3_apply_real_quadratic(vec3 k, clifford3 a, double q,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    if (hermite_quadratic_apply(a.data[0], q, 0,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return clifford3_center_slice(k, a, value.even, 0, value.odd, 0);
    }
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* I^2<0: represent the centre by complex numbers. */
static inline clifford3
clifford3_apply_complex_center(vec3 k, clifford3 a, double s, double q0, double q7,
    double (*f)(double), complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double scale = sqrt(-s);
    complex z = complex_make(a.data[0], scale*a.data[7]);
    complex q = complex_make(q0, scale*q7);
    complex root = complex_sqrt(q);
    complex even, odd;
    double radius = hypot(root.data[0], root.data[1]);
    if (radius == 0) {
        even = cf(z);
        odd = cderiv1(z);
    } else if (radius < 0.02*hermite_function_distance(z, f)) {
        hermite_complex_quadratic_near(z, root, cf, cderiv1, &even, &odd);
    } else {
        complex plus = cf(vec2_add(z, root));
        complex minus = cf(vec2_sub(z, root));
        even = vec2_add(vec2_scale(plus, 0.5), vec2_scale(minus, 0.5));
        odd = complex_div(vec2_sub(vec2_scale(plus, 0.5), vec2_scale(minus, 0.5)), root);
    }
    if (isfinite(even.data[0]) && isfinite(even.data[1]) &&
        isfinite(odd.data[0]) && isfinite(odd.data[1])) {
        return clifford3_center_slice(k, a, even.data[0], even.data[1]/scale,
            odd.data[0], odd.data[1]/scale);
    }
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* A zero projected u needs only f(w), even when f'(w) is singular. */
static inline int
clifford3_split_channel(clifford3 u, double w, double q,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), hermite_quadratic *result)
{
    if (vec8_is_real(u)) {
        result->even = f(w);
        result->odd = 0;
        return isfinite(result->even);
    }
    return hermite_quadratic_apply(w, q, 0,
        f, deriv1, deriv2, deriv3, cf, cderiv1, result);
}

/* I^2>0: evaluate the two real channels I=+/-sqrt(s), then recombine. */
static inline clifford3
clifford3_apply_split_center(vec3 k, clifford3 a, double s, double q0, double q7,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    double scale = sqrt(s);
    hermite_quadratic plus = {0}, minus = {0};
    clifford3 up = clifford3_center_slice(k, a, 0, 0, 1, 1/scale);
    clifford3 um = clifford3_center_slice(k, a, 0, 0, 1, -1/scale);
    int plus_ok = clifford3_split_channel(up, a.data[0] + scale*a.data[7], q0 + scale*q7,
        f, deriv1, deriv2, deriv3, cf, cderiv1, &plus);
    int minus_ok = clifford3_split_channel(um, a.data[0] - scale*a.data[7], q0 - scale*q7,
        f, deriv1, deriv2, deriv3, cf, cderiv1, &minus);
    if (plus_ok && minus_ok) {
        return clifford3_center_slice(k, a,
            plus.even/2 + minus.even/2, (plus.even/2 - minus.even/2)/scale,
            plus.odd/2 + minus.odd/2, (plus.odd/2 - minus.odd/2)/scale);
    }
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* I^2=0: differentiate the real quadratic coefficients in the I direction.
 * For q0=0 this is exactly the cubic Taylor formula, since (a-w)^4=0. */
static inline clifford3
clifford3_apply_dual_center(vec3 k, clifford3 a, double q0, double q7,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    hermite_quadratic value = {0};
    double t = a.data[7];
    if (hermite_quadratic_apply(a.data[0], q0, 1,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return clifford3_center_slice(k, a, value.even,
            t*value.even_dw + q7*value.odd_dw/2, value.odd,
            t*value.odd_dw + q7*value.odd_dq);
    }
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

/* With a=z+u, z=w+t*I and u^2=q0+q7*I, f(a)=A+B*u.
 * Handle the lower-dimensional cases first, then evaluate over the
 * complex, split or dual centre according to I^2=-k0*k1*k2. */
static inline clifford3
clifford3_apply(vec3 k, clifford3 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    const double *x = a.data;
    double k0 = k.data[0], k1 = k.data[1], k2 = k.data[2];
    double s = -k0*k1*k2;
    double q0 = k0*x[1]*x[1] + k1*x[2]*x[2] + k2*x[4]*x[4]
        - k0*k1*x[3]*x[3] - k0*k2*x[5]*x[5] - k1*k2*x[6]*x[6];
    double q7 = 2*(x[1]*x[6] - x[2]*x[5] + x[4]*x[3]);
    clifford3_center_coefficients value;

    if (vec8_is_real(a)) {
        return clifford3_real_slice(k, cf(complex_make(x[0], 0)));
    }
    if (x[7] == 0 && q7 == 0) {
        return clifford3_apply_real_quadratic(k, a, q0,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    if (clifford3_is_central(a)) {
        return clifford3_apply_central(a, s,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    if (clifford3_center_is_near_dual(x[0], x[7], q0, q7, s, f) &&
        clifford3_near_dual_coefficients(x[0], x[7], q0, q7, s,
            f, deriv1, deriv2, deriv3, cf, cderiv1, &value)) {
        return clifford3_center_slice(k, a,
            value.even_re, value.even_i, value.odd_re, value.odd_i);
    }
    if (s < 0) {
        return clifford3_apply_complex_center(k, a, s, q0, q7, f, cf, cderiv1);
    }
    if (s > 0) {
        return clifford3_apply_split_center(k, a, s, q0, q7,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    if (s == 0) {
        return clifford3_apply_dual_center(k, a, q0, q7,
            f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    return clifford3_make(NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN);
}

static inline clifford3
clifford3_exp(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_exp, double_exp_deriv1,
        double_exp_deriv2, double_exp_deriv3, complex_exp, complex_exp_deriv1);
}

static inline clifford3
clifford3_log(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_log, double_log_deriv1,
        double_log_deriv2, double_log_deriv3, complex_log, complex_log_deriv1);
}

static inline clifford3
clifford3_sqrt(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_sqrt, double_sqrt_deriv1,
        double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt, complex_sqrt_deriv1);
}

static inline clifford3
clifford3_cos(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_cos, double_cos_deriv1,
        double_cos_deriv2, double_cos_deriv3, complex_cos, complex_cos_deriv1);
}

static inline clifford3
clifford3_sin(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_sin, double_sin_deriv1,
        double_sin_deriv2, double_sin_deriv3, complex_sin, complex_sin_deriv1);
}

static inline clifford3
clifford3_tan(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_tan, double_tan_deriv1,
        double_tan_deriv2, double_tan_deriv3, complex_tan, complex_tan_deriv1);
}

static inline clifford3
clifford3_acos(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_acos, double_acos_deriv1,
        double_acos_deriv2, double_acos_deriv3, complex_acos, complex_acos_deriv1);
}

static inline clifford3
clifford3_asin(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_asin, double_asin_deriv1,
        double_asin_deriv2, double_asin_deriv3, complex_asin, complex_asin_deriv1);
}

static inline clifford3
clifford3_atan(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_atan, double_atan_deriv1,
        double_atan_deriv2, double_atan_deriv3, complex_atan, complex_atan_deriv1);
}

static inline clifford3
clifford3_cosh(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_cosh, double_cosh_deriv1,
        double_cosh_deriv2, double_cosh_deriv3, complex_cosh, complex_cosh_deriv1);
}

static inline clifford3
clifford3_sinh(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_sinh, double_sinh_deriv1,
        double_sinh_deriv2, double_sinh_deriv3, complex_sinh, complex_sinh_deriv1);
}

static inline clifford3
clifford3_tanh(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_tanh, double_tanh_deriv1,
        double_tanh_deriv2, double_tanh_deriv3, complex_tanh, complex_tanh_deriv1);
}

static inline clifford3
clifford3_acosh(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_acosh, double_acosh_deriv1,
        double_acosh_deriv2, double_acosh_deriv3, complex_acosh, complex_acosh_deriv1);
}

static inline clifford3
clifford3_asinh(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_asinh, double_asinh_deriv1,
        double_asinh_deriv2, double_asinh_deriv3, complex_asinh, complex_asinh_deriv1);
}

static inline clifford3
clifford3_atanh(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_atanh, double_atanh_deriv1,
        double_atanh_deriv2, double_atanh_deriv3, complex_atanh, complex_atanh_deriv1);
}

static inline clifford3
clifford3_square(vec3 k, clifford3 a)
{
    return clifford3_mul(k, a, a);
}

static inline clifford3
clifford3_rmul(vec3 k, double a, clifford3 b)
{
    return clifford3_mul(k, clifford3_from_re(a), b);
}

static inline clifford3
clifford3_lmul(vec3 k, clifford3 a, double b)
{
    return clifford3_mul(k, a, clifford3_from_re(b));
}

static inline clifford3
clifford3_rdiv(vec3 k, double a, clifford3 b)
{
    return clifford3_div(k, clifford3_from_re(a), b);
}

static inline clifford3
clifford3_ldiv(vec3 k, clifford3 a, double b)
{
    return clifford3_div(k, a, clifford3_from_re(b));
}

/* Integer powers by squaring; negative exponents invert the base once. */
static inline clifford3
clifford3_ipow(vec3 k, clifford3 a, int64_t b)
{
    clifford3 result = clifford3_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return clifford3_mul(k, a, a);
    }
    if (b == -1) {
        return clifford3_inv(k, a);
    }
    if (b < 0) {
        a = clifford3_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = clifford3_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = clifford3_mul(k, a, a);
        }
    }
    return result;
}

static inline clifford3
clifford3_lpows(vec3 k, clifford3 a, double b)
{
    if (b == 0) {
        return clifford3_one();
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return clifford3_square(k, a);
    }
    if (vec8_is_real(a)) {
        double value = pow(a.data[0], b);
        if (!isnan(value)) {
            return clifford3_from_re(value);
        }
        return clifford3_real_slice(k, complex_lpows(complex_make(a.data[0], 0), b));
    }
    return clifford3_exp(k, vec8_scale(clifford3_log(k, a), b));
}

static inline clifford3
clifford3_rpow(vec3 k, clifford3 a, clifford3 b)
{
    if (vec8_is_real(b)) {
        return clifford3_lpows(k, a, b.data[0]);
    }
    return clifford3_exp(k, clifford3_mul(k, clifford3_log(k, a), b));
}

static inline clifford3
clifford3_pow(vec3 k, clifford3 a, clifford3 b)
{
    return clifford3_rpow(k, a, b);
}

static inline clifford3
clifford3_lpow(vec3 k, clifford3 a, clifford3 b)
{
    return clifford3_exp(k, clifford3_mul(k, b, clifford3_log(k, a)));
}

static inline clifford3
clifford3_rpows(vec3 k, double a, clifford3 b)
{
    return clifford3_rpow(k, clifford3_from_re(a), b);
}

static inline clifford3
clifford3_expm1(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_expm1, double_expm1_deriv1,
        double_expm1_deriv2, double_expm1_deriv3, complex_expm1, complex_expm1_deriv1);
}

static inline clifford3
clifford3_log1p(vec3 k, clifford3 a)
{
    return clifford3_apply(k, a, double_log1p, double_log1p_deriv1,
        double_log1p_deriv2, double_log1p_deriv3, complex_log1p, complex_log1p_deriv1);
}

static inline clifford3
clifford3_mix(vec3 k, clifford3 a, clifford3 b, double t)
{
    clifford3 delta = clifford3_mul(k, clifford3_inv(k, a), b);
    return clifford3_mul(k, a, clifford3_lpows(k, delta, t));
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * Column j of lmat(a) contains a*e_j; row i of rmat(a) contains e_i*a. */
static inline void
clifford3_lmat(vec3 k, clifford3 a, double result[ALGPY_CLIFFORD3_N*ALGPY_CLIFFORD3_N])
{
    clifford3 basis = {{0}};
    clifford3 product;
    int row, col;
    for (col = 0; col < ALGPY_CLIFFORD3_N; col++) {
        basis.data[col] = 1;
        product = clifford3_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CLIFFORD3_N; row++) {
            result[row * ALGPY_CLIFFORD3_N + col] = product.data[row];
        }
    }
}

static inline void
clifford3_rmat(vec3 k, clifford3 a, double result[ALGPY_CLIFFORD3_N*ALGPY_CLIFFORD3_N])
{
    clifford3 basis = {{0}};
    clifford3 product;
    int row, col;
    for (row = 0; row < ALGPY_CLIFFORD3_N; row++) {
        basis.data[row] = 1;
        product = clifford3_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CLIFFORD3_N; col++) {
            result[row * ALGPY_CLIFFORD3_N + col] = product.data[col];
        }
    }
}


#endif /* ALGPY_CLIFFORD3_H */
