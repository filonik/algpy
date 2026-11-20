#ifndef ALGPY_HERMITE_COMPLEX_H
#define ALGPY_HERMITE_COMPLEX_H

#include "hermite.h"
#include "poly4_complex.h"

/* Complex values on real interpolation nodes, in ascending powers. */
static inline void
hermite_complex_power(int n, const double nodes[5], complex values[5], complex result[5])
{
    int i, j, degree = 0;
    for (j = 1; j < n; j++) {
        for (i = n - 1; i >= j; i--) {
            values[i] = vec2_scale(vec2_sub(values[i], values[i - 1]), 1/(nodes[i] - nodes[i - j]));
        }
    }
    result[0] = values[n - 1];
    for (i = n - 2; i >= 0; i--) {
        result[degree + 1] = result[degree];
        for (j = degree; j > 0; j--) {
            result[j] = vec2_sub(result[j - 1], vec2_scale(result[j], nodes[i]));
        }
        result[0] = vec2_sub(values[i], vec2_scale(result[0], nodes[i]));
        degree++;
    }
}

/* The same bounded derivative stencil as the real Hermite engine, now
 * centred at a complex point. Simplex moments retain every node separation. */
static inline int
hermite_complex_difference_stencil(int order, const complex *nodes, complex w,
    complex center, double radius, complex (*derivative)(complex), complex *result, double *bound)
{
    double even_nodes[5] = {0, 1, 4, 9, 16}, odd_nodes[5] = {1, 4, 9, 16, 0};
    complex even_values[5], odd_values[5], even[5], odd[5], coefficients[9];
    complex moments[9] = {{{1, 0}}};
    complex z = vec2_add(w, center), value = derivative(z);
    int i, j;
    if (!hermite_complex_is_finite(value)) return 0;
    *bound = poly_complex_abs(value);
    even_values[0] = complex_from_re(0);
    for (i = 1; i <= 4; i++) {
        complex p = derivative(complex_ladd(z, i*radius));
        complex m = derivative(complex_lsub(z, i*radius));
        if (!hermite_complex_is_finite(p) || !hermite_complex_is_finite(m)) return 0;
        *bound = fmax(*bound, fmax(poly_complex_abs(p), poly_complex_abs(m)));
        even_values[i] = vec2_sub(vec2_scale(vec2_add(p, m), .5), value);
        odd_values[i - 1] = vec2_scale(vec2_sub(p, m), .5/i);
    }
    hermite_complex_power(5, even_nodes, even_values, even);
    hermite_complex_power(4, odd_nodes, odd_values, odd);
    coefficients[0] = value;
    for (i = 1; i <= 4; i++) {
        coefficients[2*i] = even[i];
        coefficients[2*i - 1] = odd[i - 1];
    }
    for (i = 0; i <= order; i++) {
        complex u = vec2_scale(vec2_sub(nodes[i], center), 1/radius);
        for (j = 1; j <= 8; j++) moments[j] = vec2_add(moments[j], complex_mul(u, moments[j - 1]));
    }
    *result = complex_from_re(0);
    for (i = 0; i <= 8; i++) {
        double factorial = 1;
        for (j = 1; j <= order; j++) factorial *= i + j;
        *result = vec2_add(*result, vec2_scale(complex_mul(moments[i], coefficients[i]), 1/factorial));
    }
    return hermite_complex_is_finite(*result);
}

static inline int
hermite_complex_difference_near(int order, const complex *nodes, complex w,
    double (*f)(double), complex (*deriv2)(complex), complex (*deriv3)(complex), complex *result)
{
    complex center = {{0, 0}};
    double radius = 0, distance, bound;
    complex (*derivative)(complex) = order == 2 ? deriv2 : deriv3;
    for (int i = 0; i <= order; i++) center = vec2_add(center, vec2_scale(nodes[i], 1./(order + 1)));
    for (int i = 0; i <= order; i++) radius = fmax(radius, poly_complex_abs(vec2_sub(nodes[i], center)));
    distance = hermite_function_distance(vec2_add(w, center), f);
    if (!(radius > 0 && (radius < .01*distance ||
            (distance < 0 && radius < .01*fmax(1, poly_complex_abs(vec2_add(w, center))))))) return 0;
    if (!hermite_complex_difference_stencil(order, nodes, w, center, radius, derivative, result, &bound)) return 0;
    if (distance < 0) {
        complex other;
        double other_bound;
        if (!hermite_complex_difference_stencil(order, nodes, w, center, 2*radius, derivative, &other, &other_bound)) return 0;
        return poly_complex_abs(vec2_sub(*result, other)) <= 64*DBL_EPSILON*(bound + other_bound);
    }
    return 1;
}

/* Up to four arbitrary complex nodes, including multiplicities. The
 * callbacks are analytic continuations with matching branch choices.
 * f identifies a native function's cuts/poles; custom functions return the
 * usual negative distance and require checked derivative interpolation. */
static inline int
hermite_complex_interpolate(int n, complex w, const complex *input,
    double (*f)(double), complex (*cf)(complex), complex (*deriv1)(complex),
    complex (*deriv2)(complex), complex (*deriv3)(complex), complex result[4])
{
    complex nodes[4], differences[4];
    int i, j, degree = 0;
    if (n < 1 || n > 4) return 0;
    for (i = 0; i < n; i++) {
        nodes[i] = input[i];
        if (!hermite_complex_is_finite(nodes[i])) return 0;
    }
    for (i = 1; i < n; i++) {
        int nearest = i;
        double distance = INFINITY;
        for (j = i; j < n; j++) {
            double gap = poly_complex_abs(vec2_sub(nodes[j], nodes[i - 1]));
            if (gap < distance) { distance = gap; nearest = j; }
        }
        complex temporary = nodes[i];
        nodes[i] = nodes[nearest];
        nodes[nearest] = temporary;
    }
    for (i = 0; i < n; i++) {
        differences[i] = cf(vec2_add(w, nodes[i]));
        if (!hermite_complex_is_finite(differences[i])) return 0;
    }
    for (j = 1; j < n; j++) {
        for (i = n - 1; i >= j; i--) {
            complex gap = vec2_sub(nodes[i], nodes[i - j]);
            complex z = vec2_add(w, nodes[i]), candidate;
            if (poly_complex_abs(gap) == 0) {
                differences[i] = j == 1 ? deriv1(z) :
                    j == 2 ? vec2_scale(deriv2(z), .5) : vec2_scale(deriv3(z), 1./6);
            } else if (j == 1) {
                complex midpoint = vec2_add(w, vec2_scale(vec2_add(nodes[i - 1], nodes[i]), .5));
                double radius = poly_complex_abs(gap)/2;
                double distance = hermite_function_distance(midpoint, f);
                if (radius < .01*distance) {
                    differences[i] = hermite_first_difference(midpoint, vec2_scale(gap, .5), deriv1);
                } else if (distance < 0 && radius < .01*fmax(1, poly_complex_abs(midpoint)) &&
                        hermite_first_difference_checked(midpoint, vec2_scale(gap, .5), deriv1, &candidate)) {
                    differences[i] = candidate;
                } else {
                    differences[i] = complex_div(vec2_sub(differences[i], differences[i - 1]), gap);
                }
            } else if (hermite_complex_difference_near(j, nodes + i - j, w, f, deriv2, deriv3, &candidate)) {
                differences[i] = candidate;
            } else {
                differences[i] = complex_div(vec2_sub(differences[i], differences[i - 1]), gap);
            }
            if (!hermite_complex_is_finite(differences[i])) return 0;
        }
    }
    for (i = 0; i < 4; i++) result[i] = complex_from_re(0);
    result[0] = differences[n - 1];
    for (i = n - 2; i >= 0; i--) {
        result[degree + 1] = result[degree];
        for (j = degree; j > 0; j--) result[j] = vec2_sub(result[j - 1], complex_mul(nodes[i], result[j]));
        result[0] = vec2_sub(differences[i], complex_mul(nodes[i], result[0]));
        degree++;
    }
    for (i = 0; i < n; i++) if (!hermite_complex_is_finite(result[i])) return 0;
    return 1;
}

static inline int
poly4_complex_apply(const complex p[5], complex w,
    double (*f)(double), complex (*cf)(complex), complex (*deriv1)(complex),
    complex (*deriv2)(complex), complex (*deriv3)(complex), complex result[4])
{
    complex roots[4];
    poly4_complex_coeffs_to_roots(p, roots);
    return hermite_complex_interpolate(4, w, roots, f, cf, deriv1, deriv2, deriv3, result);
}

#endif /* ALGPY_HERMITE_COMPLEX_H */
