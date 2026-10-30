#ifndef ALGPY_HERMITE_H
#define ALGPY_HERMITE_H

#include <float.h>

#include "hermite_quadratic.h"
#include "poly4.h"

/* Small functional calculus for real polynomials of degree at most four.
 * Nodes are offsets from w; output coefficients are in the same variable.
 * Repeated nodes impose ordinary Hermite derivative conditions. */

static inline int
hermite_complex_is_finite(complex z)
{
    return isfinite(z.data[0]) && isfinite(z.data[1]);
}

/* Convert a real Newton interpolant to ascending power coefficients. */
static inline void
hermite_real_power(int n, const double nodes[5], double values[5], double result[5])
{
    int i, j, degree = 0;
    for (j = 1; j < n; j++) {
        for (i = n - 1; i >= j; i--) {
            values[i] = (values[i] - values[i - 1])/(nodes[i] - nodes[i - j]);
        }
    }
    result[0] = values[n - 1];
    for (i = n - 2; i >= 0; i--) {
        result[degree + 1] = result[degree];
        for (j = degree; j > 0; j--) {
            result[j] = result[j - 1] - nodes[i]*result[j];
        }
        result[0] = values[i] - nodes[i]*result[0];
        degree++;
    }
}

/* Check four-point quadrature against eight-point quadrature for custom
 * callbacks. A rapidly varying function can fail this check even on a very
 * short interval; in that case ordinary divided differences remain in use. */
static inline int
hermite_first_difference_checked(complex z, complex root,
    complex (*deriv1)(complex), complex *result)
{
    static const double nodes[4] = {.1834346424956498049, .5255324099163289858,
        .7966664774136267396, .9602898564975362317};
    static const double weights[4] = {.1813418916891809915, .1568533229389436437,
        .1111905172266872353, .0506142681451881296};
    complex coarse = hermite_first_difference(z, root, deriv1);
    double bound = 0;
    int i;
    *result = complex_from_re(0);
    for (i = 0; i < 4; i++) {
        complex offset = vec2_scale(root, nodes[i]);
        complex p = deriv1(vec2_add(z, offset)), m = deriv1(vec2_sub(z, offset));
        *result = vec2_add(*result, vec2_scale(vec2_add(p, m), weights[i]));
        bound += weights[i]*(hypot(p.data[0], p.data[1]) + hypot(m.data[0], m.data[1]));
    }
    complex error = vec2_sub(*result, coarse);
    return hermite_complex_is_finite(*result) &&
        hypot(error.data[0], error.data[1]) <= 64*DBL_EPSILON*bound;
}

/* A divided difference of order j equals the integral of f^(j) over a
 * simplex. Interpolate that derivative on the bounded
 * stencil 0,+/-r,...,+/-4r, then integrate its polynomial exactly. Complete
 * homogeneous polynomials give the simplex moments; no small gap is divided
 * out. All coordinates are scaled by r, including the simplex moments. */
static inline int
hermite_difference_stencil(int order, const complex *nodes, double w, double center, double radius,
    double (*derivative)(double), complex *result, double *bound)
{
    double value, factorial;
    double even_nodes[5] = {0, 1, 4, 9, 16};
    double odd_nodes[5] = {1, 4, 9, 16, 0};
    double even_values[5], odd_values[5], even[5], odd[5], coefficients[9];
    complex moments[9] = {{{1, 0}}};
    int i, j;
    value = derivative(w + center);
    if (!isfinite(value)) {
        return 0;
    }
    *bound = fabs(value);
    even_values[0] = 0;
    for (i = 1; i <= 4; i++) {
        double p = derivative(w + center + i*radius);
        double m = derivative(w + center - i*radius);
        if (!isfinite(p) || !isfinite(m)) {
            return 0;
        }
        *bound = fmax(*bound, fmax(fabs(p), fabs(m)));
        even_values[i] = p/2 + m/2 - value;
        odd_values[i - 1] = (p/2 - m/2)/i;
    }
    hermite_real_power(5, even_nodes, even_values, even);
    hermite_real_power(4, odd_nodes, odd_values, odd);
    coefficients[0] = value;
    for (i = 1; i <= 4; i++) {
        coefficients[2*i] = even[i];
        coefficients[2*i - 1] = odd[i - 1];
    }
    for (i = 0; i <= order; i++) {
        complex u = complex_make((nodes[i].data[0] - center)/radius, nodes[i].data[1]/radius);
        for (j = 1; j <= 8; j++) {
            moments[j] = vec2_add(moments[j], complex_mul(u, moments[j - 1]));
        }
    }
    *result = complex_from_re(0);
    for (i = 0; i <= 8; i++) {
        factorial = 1;
        for (j = 1; j <= order; j++) {
            factorial *= i + j;
        }
        *result = vec2_add(*result, vec2_scale(moments[i], coefficients[i]/factorial));
    }
    return hermite_complex_is_finite(*result);
}

/* Native functions use their known distance to a branch cut or pole.
 * Custom callbacks need agreement between two stencil radii instead.
 * A zero distance is a native cut/pole and must never use quadrature. */
static inline int
hermite_difference_near(int order, const complex *nodes, double w,
    double (*f)(double), double (*deriv2)(double), double (*deriv3)(double),
    complex *result, double *bound)
{
    double center = 0, radius = 0, distance;
    double (*derivative)(double) = order == 2 ? deriv2 : deriv3;
    int i;
    for (i = 0; i <= order; i++) {
        center += nodes[i].data[0]/(order + 1);
    }
    for (i = 0; i <= order; i++) {
        radius = fmax(radius, hypot(nodes[i].data[0] - center, nodes[i].data[1]));
    }
    distance = hermite_function_distance(complex_make(w + center, 0), f);
    if (!(radius > 0 && (radius < 0.01*distance ||
            (distance < 0 && radius < 0.01*fmax(1, fabs(w + center)))))) {
        return 0;
    }
    if (!hermite_difference_stencil(order, nodes, w, center, radius, derivative, result, bound)) {
        return 0;
    }
    if (distance < 0) {
        complex other, error;
        double other_bound;
        if (!hermite_difference_stencil(order, nodes, w, center, 2*radius, derivative, &other, &other_bound)) {
            return 0;
        }
        error = vec2_sub(*result, other);
        return hypot(error.data[0], error.data[1]) <= 64*DBL_EPSILON*(*bound + other_bound);
    }
    return 1;
}

/* Real/conjugate nodes, including multiplicities, in arbitrary order.
 * Return 0 for an unsupported domain or a nonreal interpolating polynomial.
 * Do not merge nearby nodes: their separation is retained in the moments. */
static inline int
hermite_interpolate(int n, double w, const complex *input,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), double result[4])
{
    complex nodes[4], differences[4], coefficients[4];
    double bounds[4], difference_bounds[4];
    int i, j, degree = 0;
    if (n < 1 || n > 4) {
        return 0;
    }
    for (i = 0; i < n; i++) {
        nodes[i] = input[i];
        if (!hermite_complex_is_finite(nodes[i])) {
            return 0;
        }
    }
    /* Group nearby nodes before building the table. Ordering only by real
     * part can interleave two nearly double complex conjugate pairs, placing
     * a tiny denominator around a large intervening spectral gap. */
    for (i = 1; i < n; i++) {
        int nearest = i;
        double distance = INFINITY;
        for (j = i; j < n; j++) {
            double gap = hypot(nodes[j].data[0] - nodes[i - 1].data[0],
                               nodes[j].data[1] - nodes[i - 1].data[1]);
            if (gap < distance) {
                distance = gap;
                nearest = j;
            }
        }
        complex temporary = nodes[i];
        nodes[i] = nodes[nearest];
        nodes[nearest] = temporary;
    }
    for (i = 0; i < n; i++) {
        complex z = complex_make(w + nodes[i].data[0], nodes[i].data[1]);
        differences[i] = z.data[1] == 0 ? complex_from_re(f(z.data[0])) : cf(z);
        if (!hermite_complex_is_finite(differences[i])) {
            return 0;
        }
        difference_bounds[i] = hypot(differences[i].data[0], differences[i].data[1]);
    }
    for (j = 1; j < n; j++) {
        for (i = n - 1; i >= j; i--) {
            complex gap = vec2_sub(nodes[i], nodes[i - j]);
            complex z = complex_make(w + nodes[i].data[0], nodes[i].data[1]);
            complex candidate;
            double bound;
            if (gap.data[0] == 0 && gap.data[1] == 0) {
                if (j == 1) {
                    differences[i] = z.data[1] == 0 ? complex_from_re(deriv1(z.data[0])) : cderiv1(z);
                } else if (z.data[1] == 0) {
                    differences[i] = complex_from_re(j == 2 ? deriv2(z.data[0])/2 : deriv3(z.data[0])/6);
                } else {
                    return 0; /* A real quartic has no triple nonreal root. */
                }
                bound = hypot(differences[i].data[0], differences[i].data[1]);
            } else if (j == 1) {
                complex midpoint = complex_make(w + nodes[i - 1].data[0]/2 + nodes[i].data[0]/2,
                    nodes[i - 1].data[1]/2 + nodes[i].data[1]/2);
                double radius = hypot(gap.data[0], gap.data[1])/2;
                double distance = hermite_function_distance(midpoint, f);
                if (radius < 0.01*distance) {
                    differences[i] = hermite_first_difference(midpoint, vec2_scale(gap, .5), cderiv1);
                    bound = hypot(differences[i].data[0], differences[i].data[1]);
                } else if (distance < 0 && radius < 0.01*fmax(1, hypot(midpoint.data[0], midpoint.data[1])) &&
                        hermite_first_difference_checked(midpoint, vec2_scale(gap, .5), cderiv1, &candidate)) {
                    differences[i] = candidate;
                    bound = hypot(differences[i].data[0], differences[i].data[1]);
                } else {
                    bound = (difference_bounds[i] + difference_bounds[i - 1])/(2*radius);
                    differences[i] = complex_div(vec2_sub(differences[i], differences[i - 1]), gap);
                }
            } else {
                if (hermite_difference_near(j, nodes + i - j, w, f, deriv2, deriv3, &candidate, &bound)) {
                    differences[i] = candidate;
                } else {
                    bound = (difference_bounds[i] + difference_bounds[i - 1])/hypot(gap.data[0], gap.data[1]);
                    differences[i] = complex_div(vec2_sub(differences[i], differences[i - 1]), gap);
                }
            }
            if (!hermite_complex_is_finite(differences[i])) {
                return 0;
            }
            difference_bounds[i] = bound;
        }
    }
    coefficients[0] = differences[n - 1];
    bounds[0] = difference_bounds[n - 1];
    for (i = n - 2; i >= 0; i--) {
        double radius = hypot(nodes[i].data[0], nodes[i].data[1]);
        coefficients[degree + 1] = coefficients[degree];
        bounds[degree + 1] = bounds[degree];
        for (j = degree; j > 0; j--) {
            coefficients[j] = vec2_sub(coefficients[j - 1], complex_mul(nodes[i], coefficients[j]));
            bounds[j] = bounds[j - 1] + radius*bounds[j];
        }
        coefficients[0] = vec2_sub(differences[i], complex_mul(nodes[i], coefficients[0]));
        bounds[0] = difference_bounds[i] + radius*bounds[0];
        degree++;
    }
    for (i = 0; i < 4; i++) {
        result[i] = 0;
        if (i < n) {
            if (!hermite_complex_is_finite(coefficients[i]) ||
                fabs(coefficients[i].data[1]) > 256*DBL_EPSILON*bounds[i] + 16*DBL_MIN) {
                return 0;
            }
            result[i] = coefficients[i].data[0];
        }
    }
    return 1;
}

/* Solve a real quartic, then construct the cubic Hermite interpolant for
 * f(w+x). Power-of-two scaling keeps Ferrari's intermediate coefficients
 * moderate. Scaling the variable preserves small root separations. */
static inline int
poly4_apply(const double p[5], double w,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex), double result[4])
{
    double normalized[5], radius;
    complex roots[4];
    int i, exponent = 0;
    for (i = 0; i < 5; i++) {
        normalized[i] = p[i]/p[4];
        if (!isfinite(normalized[i])) {
            return 0;
        }
    }
    radius = fmax(fabs(normalized[3]), fmax(sqrt(fabs(normalized[2])),
        fmax(cbrt(fabs(normalized[1])), sqrt(sqrt(fabs(normalized[0]))))));
    if (radius > 0) {
        frexp(radius, &exponent);
        for (i = 0; i < 4; i++) {
            normalized[i] = scalbn(normalized[i], -(4 - i)*exponent);
        }
    }
    poly4_coeffs_to_roots(normalized, roots);
    for (i = 0; i < 4; i++) {
        roots[i].data[0] = scalbn(roots[i].data[0], exponent);
        roots[i].data[1] = scalbn(roots[i].data[1], exponent);
    }
    return hermite_interpolate(4, w, roots, f, deriv1, deriv2, deriv3, cf, cderiv1, result);
}

#endif /* ALGPY_HERMITE_H */
