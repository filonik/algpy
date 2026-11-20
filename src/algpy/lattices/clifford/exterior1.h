/* Exterior algebra on 1 generator as a lattice.
 * The Clifford algebra with every generator square 0; 2 real coefficients.
 * join is the exterior product, with identity bottom = 1.
 * meet is the regressive product, with identity top = e0.
 * rdual(e_a) join e_a = top and e_a join ldual(e_a) = top for basis blades,
 * and meet(a, b) = rdual(join(ldual(a), ldual(b))).
 * In odd dimensions the two duals coincide and are also exposed as dual. */
#ifndef ALGPY_EXTERIOR1_H
#define ALGPY_EXTERIOR1_H

#include <string.h>

#include "../../vectorspaces/vec2.h"

#define ALGPY_EXTERIOR1_N 2

typedef vec2 exterior1;

static inline exterior1
exterior1_make(const double coeffs[ALGPY_EXTERIOR1_N])
{
    exterior1 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline exterior1
exterior1_grade0(exterior1 a)
{
    const double *x = a.data;
    return exterior1_make((const double[]){x[0], 0});
}

static inline exterior1
exterior1_grade1(exterior1 a)
{
    const double *x = a.data;
    return exterior1_make((const double[]){0, x[1]});
}

static inline exterior1
exterior1_from_re(double a)
{
    return exterior1_make((const double[ALGPY_EXTERIOR1_N]){a});
}

static inline double
exterior1_to_re(exterior1 a)
{
    return a.data[0];
}

static inline exterior1
exterior1_zero(void)
{
    return vec2_zero();
}

static inline exterior1
exterior1_bottom(void)
{
    return exterior1_make((const double[ALGPY_EXTERIOR1_N]){[0] = 1});
}

static inline exterior1
exterior1_top(void)
{
    return exterior1_make((const double[ALGPY_EXTERIOR1_N]){[1] = 1});
}

static inline exterior1
exterior1_pos(exterior1 a)
{
    return vec2_pos(a);
}

static inline exterior1
exterior1_neg(exterior1 a)
{
    return vec2_neg(a);
}

static inline exterior1
exterior1_add(exterior1 a, exterior1 b)
{
    return vec2_add(a, b);
}

static inline exterior1
exterior1_sub(exterior1 a, exterior1 b)
{
    return vec2_sub(a, b);
}

static inline exterior1
exterior1_rmul(double a, exterior1 b)
{
    return vec2_scale(b, a);
}

static inline exterior1
exterior1_lmul(exterior1 a, double b)
{
    return vec2_scale(a, b);
}

static inline exterior1
exterior1_ldiv(exterior1 a, double b)
{
    return vec2_div_scalar(a, b);
}

static inline int
exterior1_equal(exterior1 a, exterior1 b)
{
    return vec2_all_eq(a, b);
}

static inline int
exterior1_not_equal(exterior1 a, exterior1 b)
{
    return !exterior1_equal(a, b);
}

static inline int
exterior1_is_real(exterior1 a)
{
    return vec2_is_real(a);
}

static inline exterior1
exterior1_join(exterior1 a, exterior1 b)
{
    const double *x = a.data, *y = b.data;
    return exterior1_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0]});
}

static inline exterior1
exterior1_meet(exterior1 a, exterior1 b)
{
    const double *x = a.data, *y = b.data;
    return exterior1_make((const double[]){
        x[0]*y[1] + x[1]*y[0],
        x[1]*y[1]});
}

static inline exterior1
exterior1_rdual(exterior1 a)
{
    const double *x = a.data;
    return exterior1_make((const double[]){
        x[1],
        x[0]});
}

static inline exterior1
exterior1_ldual(exterior1 a)
{
    const double *x = a.data;
    return exterior1_make((const double[]){
        x[1],
        x[0]});
}

static inline exterior1
exterior1_dual(exterior1 a)
{
    return exterior1_rdual(a);
}

static inline double
exterior1_norm2(exterior1 a)
{
    double result = 0;
    int i;
    for (i = 0; i < ALGPY_EXTERIOR1_N; i++) {
        result += a.data[i]*a.data[i];
    }
    return result;
}

static inline void
exterior1_basis(exterior1 result[2])
{
    int i, j;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
exterior1_gens(exterior1 result[1])
{
    result[0] = exterior1_make((const double[ALGPY_EXTERIOR1_N]){[1] = 1});
}

#endif /* ALGPY_EXTERIOR1_H */
