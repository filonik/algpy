#ifndef ALGPY_POLY4_COMPLEX_H
#define ALGPY_POLY4_COMPLEX_H

#include "poly4.h"

static inline double
poly_complex_abs(complex z)
{
    return hypot(z.data[0], z.data[1]);
}

/* Ascending complex coefficients. These small closed-form solvers support
 * quartic functional calculus over a complex centre. */
static inline void
poly2_complex_coeffs_to_roots(const complex p[3], complex result[2])
{
    complex b = complex_div(p[1], p[2]), c = complex_div(p[0], p[2]);
    complex root = complex_sqrt(vec2_sub(complex_square(b), vec2_scale(c, 4)));
    complex plus = vec2_add(b, root), minus = vec2_sub(b, root);
    complex q = vec2_scale(poly_complex_abs(plus) >= poly_complex_abs(minus) ? plus : minus, -.5);
    result[0] = q;
    result[1] = poly_complex_abs(q) == 0 ? complex_from_re(0) : complex_div(c, q);
}

static inline void
poly3_complex_coeffs_to_roots(const complex p[4], complex result[3])
{
    complex a = complex_div(p[2], p[3]), b = complex_div(p[1], p[3]);
    complex c = complex_div(p[0], p[3]), shift = vec2_scale(a, 1./3);
    complex u = vec2_scale(vec2_sub(b, complex_mul(a, shift)), 1./3);
    complex v = vec2_scale(vec2_add(complex_mul(vec2_sub(vec2_scale(complex_square(shift), 2), b), shift), c), .5);
    complex d = vec2_add(complex_square(v), complex_mul(complex_square(u), u));
    double scale = poly_complex_abs(v)*poly_complex_abs(v) + pow(poly_complex_abs(u), 3);
    if (poly_complex_abs(c) == 0) {
        complex quadratic[3] = {b, a, {{1, 0}}};
        result[0] = complex_from_re(0);
        poly2_complex_coeffs_to_roots(quadratic, result + 1);
    } else if (poly_complex_abs(d) <= 8*DBL_EPSILON*scale) {
        complex repeated = poly_complex_abs(u) == 0 ? complex_from_re(0) : complex_div(v, u);
        result[0] = vec2_sub(vec2_scale(repeated, -2), shift);
        result[1] = result[2] = vec2_sub(repeated, shift);
    } else {
        complex root = complex_sqrt(d);
        complex plus = vec2_add(vec2_neg(v), root), minus = vec2_sub(vec2_neg(v), root);
        complex r = complex_exp(vec2_scale(complex_log(
            poly_complex_abs(plus) >= poly_complex_abs(minus) ? plus : minus), 1./3));
        complex s = poly_complex_abs(r) == 0 ? complex_from_re(0) : complex_div(vec2_neg(u), r);
        complex omega = complex_make((const double[]){-.5, sqrt(3.)/2});
        complex other = complex_conjugate(omega);
        result[0] = vec2_sub(vec2_add(r, s), shift);
        result[1] = vec2_sub(vec2_add(complex_mul(omega, r), complex_mul(other, s)), shift);
        result[2] = vec2_sub(vec2_add(complex_mul(other, r), complex_mul(omega, s)), shift);
    }
}

/* Power-of-two scaling bounds Ferrari's intermediates. Real polynomials
 * retain the existing real solver, including its repeated-root cases. */
static inline void
poly4_complex_coeffs_to_roots(const complex p[5], complex result[4])
{
    complex normalized[5];
    double radius = 0;
    int exponent = 0, real = 1, i;
    for (i = 0; i < 5; i++) {
        normalized[i] = complex_div(p[i], p[4]);
        if (!isfinite(normalized[i].data[0]) || !isfinite(normalized[i].data[1])) {
            for (int j = 0; j < 4; j++) result[j] = complex_make((const double[]){NAN, NAN});
            return;
        }
        if (normalized[i].data[1] != 0) real = 0;
    }
    radius = fmax(poly_complex_abs(normalized[3]), fmax(sqrt(poly_complex_abs(normalized[2])),
        fmax(cbrt(poly_complex_abs(normalized[1])), sqrt(sqrt(poly_complex_abs(normalized[0]))))));
    if (radius > 0) frexp(radius, &exponent);
    for (i = 0; i < 4; i++) {
        normalized[i].data[0] = scalbn(normalized[i].data[0], -(4 - i)*exponent);
        normalized[i].data[1] = scalbn(normalized[i].data[1], -(4 - i)*exponent);
    }
    if (real) {
        double coefficients[5];
        for (i = 0; i < 5; i++) coefficients[i] = normalized[i].data[0];
        poly4_coeffs_to_roots(coefficients, result);
    } else {
        complex a = normalized[3], b = normalized[2], c = normalized[1], d = normalized[0];
        complex shift = vec2_scale(a, .25), shift2 = complex_square(shift);
        complex pp = vec2_sub(b, vec2_scale(shift2, 6));
        complex qq = vec2_add(complex_mul(vec2_sub(vec2_scale(shift2, 8), vec2_scale(b, 2)), shift), c);
        complex rr = vec2_add(complex_mul(vec2_sub(complex_mul(vec2_sub(b, vec2_scale(shift2, 3)), shift), c), shift), d);
        double qscale = 8*pow(poly_complex_abs(shift), 3) + 2*poly_complex_abs(b)*poly_complex_abs(shift) + poly_complex_abs(c);
        qscale += pow(poly_complex_abs(pp), 1.5) + pow(poly_complex_abs(rr), .75);
        if (poly_complex_abs(d) == 0) {
            poly3_complex_coeffs_to_roots(normalized + 1, result);
            result[3] = complex_from_re(0);
        } else if (poly_complex_abs(qq) <= 4*DBL_EPSILON*qscale) {
            complex quadratic[3] = {rr, pp, {{1, 0}}}, squares[2];
            poly2_complex_coeffs_to_roots(quadratic, squares);
            complex u = complex_sqrt(squares[0]), v = complex_sqrt(squares[1]);
            result[0] = vec2_sub(u, shift);
            result[1] = vec2_sub(vec2_neg(u), shift);
            result[2] = vec2_sub(v, shift);
            result[3] = vec2_sub(vec2_neg(v), shift);
        } else {
            complex cubic[4] = {vec2_neg(complex_square(qq)),
                vec2_sub(complex_square(pp), vec2_scale(rr, 4)), vec2_scale(pp, 2), {{1, 0}}};
            complex roots[3], m, s, ratio;
            poly3_complex_coeffs_to_roots(cubic, roots);
            m = roots[0];
            for (i = 1; i < 3; i++) if (poly_complex_abs(roots[i]) > poly_complex_abs(m)) m = roots[i];
            s = complex_sqrt(m);
            ratio = complex_div(qq, s);
            complex left[3] = {vec2_scale(vec2_sub(vec2_add(pp, m), ratio), .5), s, {{1, 0}}};
            complex right[3] = {vec2_scale(vec2_add(vec2_add(pp, m), ratio), .5), vec2_neg(s), {{1, 0}}};
            poly2_complex_coeffs_to_roots(left, result);
            poly2_complex_coeffs_to_roots(right, result + 2);
            for (i = 0; i < 4; i++) result[i] = vec2_sub(result[i], shift);
        }
    }
    for (i = 0; i < 4; i++) {
        result[i].data[0] = scalbn(result[i].data[0], exponent);
        result[i].data[1] = scalbn(result[i].data[1], exponent);
    }
}

#endif /* ALGPY_POLY4_COMPLEX_H */
