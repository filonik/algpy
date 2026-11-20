#ifndef ALGPY_HERMITE8_H
#define ALGPY_HERMITE8_H

#include "hermite_complex.h"

/* Real/conjugate nodes with real multiplicity up to eight and nonreal
 * multiplicity up to four. Clifford5 and higher PGA specializations use this
 * larger interpolant; their nodes still come from quartic solvers. */

/* Cauchy coefficients of f''' on a circle inside a native analytic domain.
 * Pair conjugate samples so the coefficients remain exactly real. Using
 * f''' retains small derivatives when the value has rounded to a constant. */
static inline int
hermite8_cauchy_jets(double w, const complex *nodes, int n,
    double (*f)(double), complex (*deriv3)(complex),
    double *center, double *radius, double jets[29])
{
    double spread = 0;
    *center = 0;
    for (int i = 0; i < n; i++) *center += nodes[i].data[0]/n;
    for (int i = 0; i < n; i++) spread = fmax(spread,
        hypot(nodes[i].data[0] - *center, nodes[i].data[1]));
    double distance = hermite_function_distance(complex_make((const double[]){w + *center, 0}), f);
    /* Far along the real axis, tanh's poles are distant but its small
     * derivatives still vary on a unit exponential scale. */
    if (f == double_tanh) distance = fmin(distance, 1);
    if (!(spread > 0 && spread < .1*distance && isfinite(distance))) return 0;
    *radius = distance/2;
    for (int i = 0; i < 29; i++) jets[i] = 0;
    for (int sample = 0; sample < 32; sample++) {
        double angle = M_PI*(sample + .5)/32;
        complex unit = complex_make((const double[]){cos(angle), sin(angle)});
        complex value = deriv3(complex_make((const double[]){w + *center + *radius*unit.data[0],
                                            *radius*unit.data[1]}));
        if (!hermite_complex_is_finite(value)) return 0;
        complex phase = {{1, 0}};
        for (int i = 0; i < 29; i++) {
            jets[i] += complex_mul(value, phase).data[0]/32;
            phase = complex_mul(phase, complex_conjugate(unit));
        }
    }
    return 1;
}

/* f[x0,...,xj] = sum_m f^(j+m)(c)/(j+m)! * h_m(x-c).
 * Complete homogeneous moments are scaled by the contour radius; the
 * factorial ratio converts the sampled third derivative to f's series. */
static inline int
hermite8_cauchy_difference(int order, const complex *nodes,
    double center, double radius, const double jets[29], complex *result, double *bound)
{
    complex moments[25] = {{{1, 0}}};
    double moment_bounds[25] = {1};
    for (int i = 0; i <= order; i++) {
        complex u = complex_make((const double[]){(nodes[i].data[0] - center)/radius, nodes[i].data[1]/radius});
        double size = poly_complex_abs(u);
        for (int j = 1; j <= 24; j++) {
            moments[j] = vec2_add(moments[j], complex_mul(u, moments[j - 1]));
            moment_bounds[j] += size*moment_bounds[j - 1];
        }
    }
    *result = complex_from_re(0);
    *bound = 0;
    for (int i = 0; i <= 24; i++) {
        double degree = order + i;
        double coefficient = jets[order + i - 3]/(degree*(degree - 1)*(degree - 2));
        *result = vec2_add(*result, vec2_scale(moments[i], coefficient));
        /* Conjugate moments can cancel almost completely. Bound their
         * arithmetic before cancellation, including for odd coefficients. */
        *bound += fabs(coefficient)*moment_bounds[i];
    }
    for (int i = 3; i < order; i++) {
        *result = vec2_div_scalar(*result, radius);
        *bound /= radius;
    }
    return hermite_complex_is_finite(*result);
}

static inline int
hermite8_difference_near(int order, const complex *nodes, double w,
    double (*f)(double), double (*derivatives[7])(double),
    complex *result, double *bound)
{
    double center = 0, radius = 0, distance;
    double (*derivative)(double) = derivatives[order - 1];
    int i;
    for (i = 0; i <= order; i++) {
        center += nodes[i].data[0]/(order + 1);
    }
    for (i = 0; i <= order; i++) {
        radius = fmax(radius, hypot(nodes[i].data[0] - center, nodes[i].data[1]));
    }
    distance = hermite_function_distance(complex_make((const double[]){w + center, 0}), f);
    if (!(radius > 0 && (radius < .01*distance ||
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

static inline int
hermite8_interpolate(int n, double w, const complex *input,
    double (*f)(double), double (*derivatives[7])(double),
    complex (*cf)(complex), complex (*cderivatives[3])(complex), double result[8])
{
    complex nodes[8], differences[8], coefficients[8];
    double bounds[8], difference_bounds[8];
    double jet_center = 0, jet_radius = 1, jets[29];
    int i, j, degree = 0;
    if (n < 1 || n > 8) {
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
        complex z = complex_make((const double[]){w + nodes[i].data[0], nodes[i].data[1]});
        differences[i] = z.data[1] == 0 ? complex_from_re(f(z.data[0])) : cf(z);
        if (!hermite_complex_is_finite(differences[i])) {
            return 0;
        }
        difference_bounds[i] = hypot(differences[i].data[0], differences[i].data[1]);
    }
    int with_jets = n >= 5 && hermite8_cauchy_jets(w, nodes, n, f, cderivatives[2],
                                                 &jet_center, &jet_radius, jets);
    for (j = 1; j < n; j++) {
        for (i = n - 1; i >= j; i--) {
            complex gap = vec2_sub(nodes[i], nodes[i - j]);
            complex z = complex_make((const double[]){w + nodes[i].data[0], nodes[i].data[1]});
            complex candidate;
            double bound;
            if (gap.data[0] == 0 && gap.data[1] == 0) {
                double factorial = 1;
                for (int order = 2; order <= j; order++) factorial *= order;
                if (z.data[1] == 0) {
                    differences[i] = complex_from_re(derivatives[j - 1](z.data[0])/factorial);
                } else if (j <= 3) {
                    differences[i] = vec2_scale(cderivatives[j - 1](z), 1/factorial);
                } else {
                    return 0; /* Conjugate pairs have multiplicity at most four. */
                }
                bound = hypot(differences[i].data[0], differences[i].data[1]);
            } else if (j == 1) {
                complex midpoint = complex_make((const double[]){w + nodes[i - 1].data[0]/2 + nodes[i].data[0]/2,
                    nodes[i - 1].data[1]/2 + nodes[i].data[1]/2});
                double radius = hypot(gap.data[0], gap.data[1])/2;
                double distance = hermite_function_distance(midpoint, f);
                if (radius < 0.01*distance) {
                    differences[i] = hermite_first_difference(midpoint, vec2_scale(gap, .5), cderivatives[0]);
                    bound = hypot(differences[i].data[0], differences[i].data[1]);
                } else if (distance < 0 && radius < 0.01*fmax(1, hypot(midpoint.data[0], midpoint.data[1])) &&
                        hermite_first_difference_checked(midpoint, vec2_scale(gap, .5), cderivatives[0], &candidate)) {
                    differences[i] = candidate;
                    bound = hypot(differences[i].data[0], differences[i].data[1]);
                } else {
                    bound = (difference_bounds[i] + difference_bounds[i - 1])/(2*radius);
                    differences[i] = complex_div(vec2_sub(differences[i], differences[i - 1]), gap);
                }
            } else {
                /* A nonreal repeated root can be split by root-solver
                 * roundoff. Centre its order-2/3 stencil in the complex
                 * plane, rather than measuring its gap from the real axis. */
                if (j >= 4 && with_jets && hermite8_cauchy_difference(j, nodes + i - j,
                        jet_center, jet_radius, jets, &candidate, &bound)) {
                    differences[i] = candidate;
                } else if (j <= 3 && hermite_complex_difference_near(j, nodes + i - j,
                        complex_from_re(w), f, cderivatives[1], cderivatives[2], &candidate)) {
                    differences[i] = candidate;
                    bound = poly_complex_abs(candidate);
                } else if (hermite8_difference_near(j, nodes + i - j, w, f, derivatives, &candidate, &bound)) {
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
    for (i = 0; i < 8; i++) {
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

#endif /* ALGPY_HERMITE8_H */
