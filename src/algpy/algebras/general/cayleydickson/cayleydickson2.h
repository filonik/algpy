#ifndef ALGPY_CAYLEYDICKSON2_H
#define ALGPY_CAYLEYDICKSON2_H

#define ALGPY_CAYLEYDICKSON2_N 4

#include <math.h>

#include "../../../fields/complex.h"
#include "../../../vectorspaces/vec2.h"
#include "../../../vectorspaces/vec4.h"

typedef vec4 cayleydickson2;

static inline cayleydickson2
cayleydickson2_make(double w, double x, double y, double z)
{
    cayleydickson2 r = {{w, x, y, z}};
    return r;
}

static inline cayleydickson2
cayleydickson2_re(cayleydickson2 a)
{
    return cayleydickson2_make(a.data[0], 0, 0, 0);
}

static inline cayleydickson2
cayleydickson2_im(cayleydickson2 a)
{
    return cayleydickson2_make(0, a.data[1], a.data[2], a.data[3]);
}

static inline cayleydickson2
cayleydickson2_from_re(double a)
{
    return cayleydickson2_make(a, 0, 0, 0);
}

static inline double
cayleydickson2_to_re(cayleydickson2 a)
{
    return a.data[0];
}

static inline cayleydickson2
cayleydickson2_one(void)
{
    return cayleydickson2_from_re(1);
}

static inline cayleydickson2
cayleydickson2_mul(vec2 h, cayleydickson2 a, cayleydickson2 b)
{
    const double *k = h.data, *x = a.data, *y = b.data;
    return cayleydickson2_make(
      x[0]*y[0] + k[0]*x[1]*y[1] + k[1]*x[2]*y[2] - k[0]*k[1]*x[3]*y[3],
      x[0]*y[1] + x[1]*y[0] - k[1]*x[2]*y[3] + k[1]*x[3]*y[2],
      x[0]*y[2] + k[0]*x[1]*y[3] + x[2]*y[0] - k[0]*x[3]*y[1],
      x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0]);
}

static inline void
cayleydickson2_basis(cayleydickson2 result[ALGPY_CAYLEYDICKSON2_N])
{
    int i, j;
    for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON2_N; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
cayleydickson2_coeffs(
    vec2 h,
    cayleydickson2 result[ALGPY_CAYLEYDICKSON2_N * ALGPY_CAYLEYDICKSON2_N])
{
    cayleydickson2 basis[ALGPY_CAYLEYDICKSON2_N];
    int i, j;
    cayleydickson2_basis(basis);
    for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
        for (j = 0; j < ALGPY_CAYLEYDICKSON2_N; j++) {
            result[i * ALGPY_CAYLEYDICKSON2_N + j] =
                cayleydickson2_mul(h, basis[i], basis[j]);
        }
    }
}

static inline cayleydickson2
cayleydickson2_conjugate(cayleydickson2 a)
{
    return vec4_sub(cayleydickson2_re(a), cayleydickson2_im(a));
}

/* Real norm a*conj(a) = w^2 - k0*x^2 - k1*y^2 + k0*k1*z^2 */
static inline double
cayleydickson2_norm2(vec2 k, cayleydickson2 a)
{
    const double *x = a.data;
    return x[0]*x[0] - k.data[0]*x[1]*x[1] - k.data[1]*x[2]*x[2] +
        k.data[0]*k.data[1]*x[3]*x[3];
}

static inline cayleydickson2
cayleydickson2_inv(vec2 k, cayleydickson2 a)
{
    return vec4_div_scalar(cayleydickson2_conjugate(a), cayleydickson2_norm2(k, a));
}

static inline cayleydickson2
cayleydickson2_div(vec2 k, cayleydickson2 a, cayleydickson2 b)
{
    return cayleydickson2_mul(k, a, cayleydickson2_inv(k, b));
}

/* Write f(w + v) = real + scale*v, where v*v is a real scalar.
 * All analytic functions use the same algebra parameters as mul/inv.
 * log and sqrt use the principal branch in the subalgebra spanned by 1,v.
 * Outside that real subalgebra's domain they return four NaNs. Negative
 * real inputs use the first basis direction with negative square, if any.
 * log(0) is undefined and also returns four NaNs.
 */
static inline cayleydickson2
cayleydickson2_slice(cayleydickson2 a, double real, double scale)
{
    return cayleydickson2_make(real,
        a.data[1] == 0 ? 0 : scale*a.data[1],
        a.data[2] == 0 ? 0 : scale*a.data[2],
        a.data[3] == 0 ? 0 : scale*a.data[3]);
}

static inline cayleydickson2
cayleydickson2_negative_real(vec2 k, double real, double imag)
{
    double squares[3] = {k.data[0], k.data[1], -k.data[0]*k.data[1]};
    int i;
    for (i = 0; i < 3; ++i) {
        if (squares[i] < 0) {
            cayleydickson2 r = cayleydickson2_from_re(real);
            r.data[i + 1] = imag/sqrt(-squares[i]);
            return r;
        }
    }
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_exp(vec2 k, cayleydickson2 a)
{
    double q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    double r = sqrt(fabs(q)), e = exp(a.data[0]);
    if (q < 0) {
        return cayleydickson2_slice(a, e*cos(r), e*(sin(r)/r));
    }
    if (q > 0) {
        return cayleydickson2_slice(a, e*cosh(r), e*(sinh(r)/r));
    }
    return cayleydickson2_slice(a, e, e);
}

static inline cayleydickson2
cayleydickson2_log(vec2 k, cayleydickson2 a)
{
    double w = a.data[0];
    double q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    double r = sqrt(fabs(q));
    if (q < 0) {
        return cayleydickson2_slice(a, log(hypot(w, r)), atan2(r, w)/r);
    }
    if (q > 0 && w > r) {
        double lp = log(w + r), lm = log(w - r);
        /* log1p avoids cancellation when r is small relative to w. */
        return cayleydickson2_slice(a, 0.5*lp + 0.5*lm, (r < 0.5*w ? log1p(2*r/(w-r)) : lp-lm)/(2*r));
    }
    if (q == 0 && w > 0) {
        return cayleydickson2_slice(a, log(w), 1/w);
    }
    if (w < 0 && a.data[1] == 0 && a.data[2] == 0 && a.data[3] == 0) {
        return cayleydickson2_negative_real(k, log(-w), acos(-1.0));
    }
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_pow(vec2 k, cayleydickson2 a, cayleydickson2 b)
{
    return cayleydickson2_exp(k, cayleydickson2_mul(k, cayleydickson2_log(k, a), b));
}

static inline cayleydickson2
cayleydickson2_rpow(vec2 k, double a, cayleydickson2 b)
{
    return cayleydickson2_pow(k, cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_lpow(vec2 k, cayleydickson2 a, double b)
{
    return cayleydickson2_pow(k, a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_mix(vec2 k, cayleydickson2 a, cayleydickson2 b, double t)
{
    cayleydickson2 delta = cayleydickson2_mul(k, cayleydickson2_inv(k, a), b);
    return cayleydickson2_mul(k, a, cayleydickson2_lpow(k, delta, t));
}

static inline cayleydickson2
cayleydickson2_sqrt(vec2 k, cayleydickson2 a)
{
    double w = a.data[0];
    double q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    double r = sqrt(fabs(q));
    if (q < 0) {
        double m = hypot(w, r);
        if (w >= 0) {
            double t = sqrt(0.5*m + 0.5*w);
            return cayleydickson2_slice(a, t, 0.5/t);
        } else {
            double t = sqrt(0.5*m - 0.5*w);
            return cayleydickson2_slice(a, 0.5*(r/t), t/r);
        }
    }
    if (q > 0 && w >= r) {
        double t = 0.5*sqrt(w + r) + 0.5*sqrt(w - r);
        return cayleydickson2_slice(a, t, 0.5/t);
    }
    if (q == 0 && w > 0) {
        return cayleydickson2_slice(a, sqrt(w), 0.5/sqrt(w));
    }
    if (a.data[1] == 0 && a.data[2] == 0 && a.data[3] == 0) {
        if (w == 0)
            return cayleydickson2_from_re(0);
        if (w < 0)
            return cayleydickson2_negative_real(k, 0, sqrt(-w));
    }
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_cos(vec2 k, cayleydickson2 a)
{
    double q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    double r = sqrt(fabs(q)), w = a.data[0];
    if (q < 0) {
        return cayleydickson2_slice(a, cos(w)*cosh(r), -sin(w)*(sinh(r)/r));
    }
    if (q > 0) {
        return cayleydickson2_slice(a, cos(w)*cos(r), -sin(w)*(sin(r)/r));
    }
    return cayleydickson2_slice(a, cos(w), -sin(w));
}

static inline cayleydickson2
cayleydickson2_sin(vec2 k, cayleydickson2 a)
{
    double q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    double r = sqrt(fabs(q)), w = a.data[0];
    if (q < 0) {
        return cayleydickson2_slice(a, sin(w)*cosh(r), cos(w)*(sinh(r)/r));
    }
    if (q > 0) {
        return cayleydickson2_slice(a, sin(w)*cos(r), cos(w)*(sin(r)/r));
    }
    return cayleydickson2_slice(a, sin(w), cos(w));
}

static inline cayleydickson2
cayleydickson2_tan(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_div(k, cayleydickson2_sin(k, a), cayleydickson2_cos(k, a));
}

static inline cayleydickson2
cayleydickson2_cosh(vec2 k, cayleydickson2 a)
{
    double q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    double r = sqrt(fabs(q)), w = a.data[0];
    if (q < 0) {
        return cayleydickson2_slice(a, cosh(w)*cos(r), sinh(w)*(sin(r)/r));
    }
    if (q > 0) {
        return cayleydickson2_slice(a, cosh(w)*cosh(r), sinh(w)*(sinh(r)/r));
    }
    return cayleydickson2_slice(a, cosh(w), sinh(w));
}

static inline cayleydickson2
cayleydickson2_sinh(vec2 k, cayleydickson2 a)
{
    double q = -cayleydickson2_norm2(k, cayleydickson2_im(a));
    double r = sqrt(fabs(q)), w = a.data[0];
    if (q < 0) {
        return cayleydickson2_slice(a, sinh(w)*cos(r), cosh(w)*(sin(r)/r));
    }
    if (q > 0) {
        return cayleydickson2_slice(a, sinh(w)*cosh(r), cosh(w)*(sinh(r)/r));
    }
    return cayleydickson2_slice(a, sinh(w), cosh(w));
}

static inline cayleydickson2
cayleydickson2_tanh(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_div(k, cayleydickson2_sinh(k, a), cayleydickson2_cosh(k, a));
}

/* Minimal polynomial t^2 - 2*w*t + norm2(a) (coefficients in ascending order) */
static inline void
cayleydickson2_minpoly(vec2 k, cayleydickson2 a, double result[3])
{
    result[0] = cayleydickson2_norm2(k, a);
    result[1] = -2 * a.data[0];
    result[2] = 1;
}

/* Roots w +/- sqrt(-norm2(im(a))) of the minimal polynomial, as complex numbers. */
static inline void
cayleydickson2_eivals(vec2 k, cayleydickson2 a, complex result[2])
{
    complex s = complex_sqrt_real(-cayleydickson2_norm2(k, cayleydickson2_im(a)));
    result[0] = complex_make(a.data[0] - s.data[0], -s.data[1]);
    result[1] = complex_make(a.data[0] + s.data[0], +s.data[1]);
}

/* Common algebra API. Unfinished operations return NaN components. */
static inline cayleydickson2
cayleydickson2_zero(void)
{
    return vec4_zero();
}

static inline int
cayleydickson2_is_real(cayleydickson2 a)
{
    return vec4_is_real(a);
}

static inline cayleydickson2
cayleydickson2_pos(cayleydickson2 a)
{
    return vec4_pos(a);
}

static inline cayleydickson2
cayleydickson2_neg(cayleydickson2 a)
{
    return vec4_neg(a);
}

static inline cayleydickson2
cayleydickson2_add(cayleydickson2 a, cayleydickson2 b)
{
    return vec4_add(a, b);
}

static inline cayleydickson2
cayleydickson2_sub(cayleydickson2 a, cayleydickson2 b)
{
    return vec4_sub(a, b);
}

static inline int
cayleydickson2_equal(cayleydickson2 a, cayleydickson2 b)
{
    return vec4_all_eq(a, b);
}

static inline int
cayleydickson2_not_equal(cayleydickson2 a, cayleydickson2 b)
{
    return !cayleydickson2_equal(a, b);
}

static inline cayleydickson2
cayleydickson2_acos(vec2 k, cayleydickson2 a)
{
    /* TODO: implement cayleydickson2_acos. */
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_asin(vec2 k, cayleydickson2 a)
{
    /* TODO: implement cayleydickson2_asin. */
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_atan(vec2 k, cayleydickson2 a)
{
    /* TODO: implement cayleydickson2_atan. */
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_acosh(vec2 k, cayleydickson2 a)
{
    /* TODO: implement cayleydickson2_acosh. */
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_asinh(vec2 k, cayleydickson2 a)
{
    /* TODO: implement cayleydickson2_asinh. */
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_atanh(vec2 k, cayleydickson2 a)
{
    /* TODO: implement cayleydickson2_atanh. */
    return cayleydickson2_make(NAN, NAN, NAN, NAN);
}

static inline cayleydickson2
cayleydickson2_square(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_mul(k, a, a);
}

static inline cayleydickson2
cayleydickson2_radd(double a, cayleydickson2 b)
{
    return cayleydickson2_add(cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_ladd(cayleydickson2 a, double b)
{
    return cayleydickson2_add(a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_rsub(double a, cayleydickson2 b)
{
    return cayleydickson2_sub(cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_lsub(cayleydickson2 a, double b)
{
    return cayleydickson2_sub(a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_rmul(vec2 k, double a, cayleydickson2 b)
{
    return cayleydickson2_mul(k, cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_lmul(vec2 k, cayleydickson2 a, double b)
{
    return cayleydickson2_mul(k, a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_rdiv(vec2 k, double a, cayleydickson2 b)
{
    return cayleydickson2_div(k, cayleydickson2_from_re(a), b);
}

static inline cayleydickson2
cayleydickson2_ldiv(vec2 k, cayleydickson2 a, double b)
{
    return cayleydickson2_div(k, a, cayleydickson2_from_re(b));
}

static inline cayleydickson2
cayleydickson2_expm1(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_lsub(cayleydickson2_exp(k, a), 1);
}

static inline cayleydickson2
cayleydickson2_log1p(vec2 k, cayleydickson2 a)
{
    return cayleydickson2_log(k, cayleydickson2_ladd(a, 1));
}

static inline complex
cayleydickson2_norm(vec2 k, cayleydickson2 a)
{
    return complex_root_real(cayleydickson2_norm2(k, a), 2);
}

/* Row-major multiplication matrices: lmat acts on columns, rmat on rows.
 * With C[i,j,k] the k-th coefficient of e_i*e_j:
 * L(a)[k,j] = sum_i a[i]*C[i,j,k], R(a)[i,k] = sum_j a[j]*C[i,j,k]. */
static inline void
cayleydickson2_lmat(vec2 k, cayleydickson2 a, double result[ALGPY_CAYLEYDICKSON2_N*ALGPY_CAYLEYDICKSON2_N])
{
    cayleydickson2 coeffs[ALGPY_CAYLEYDICKSON2_N * ALGPY_CAYLEYDICKSON2_N];
    int row, col, i;
    cayleydickson2_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CAYLEYDICKSON2_N; row++) {
        for (col = 0; col < ALGPY_CAYLEYDICKSON2_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
                value += a.data[i] * coeffs[i * ALGPY_CAYLEYDICKSON2_N + col].data[row];
            }
            result[row * ALGPY_CAYLEYDICKSON2_N + col] = value;
        }
    }
}

static inline void
cayleydickson2_rmat(vec2 k, cayleydickson2 a, double result[ALGPY_CAYLEYDICKSON2_N*ALGPY_CAYLEYDICKSON2_N])
{
    cayleydickson2 coeffs[ALGPY_CAYLEYDICKSON2_N * ALGPY_CAYLEYDICKSON2_N];
    int row, col, i;
    cayleydickson2_coeffs(k, coeffs);
    for (row = 0; row < ALGPY_CAYLEYDICKSON2_N; row++) {
        for (col = 0; col < ALGPY_CAYLEYDICKSON2_N; col++) {
            double value = 0;
            for (i = 0; i < ALGPY_CAYLEYDICKSON2_N; i++) {
                value += a.data[i] * coeffs[row * ALGPY_CAYLEYDICKSON2_N + i].data[col];
            }
            result[row * ALGPY_CAYLEYDICKSON2_N + col] = value;
        }
    }
}


#endif /* ALGPY_CAYLEYDICKSON2_H */
