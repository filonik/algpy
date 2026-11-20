/* Exterior algebra on 2 generators as a lattice.
 * The Clifford algebra with every generator square 0; 4 real coefficients.
 * join is the exterior product, with identity bottom = 1.
 * meet is the regressive product, with identity top = e01.
 * rdual(e_a) join e_a = top and e_a join ldual(e_a) = top for basis blades,
 * and meet(a, b) = rdual(join(ldual(a), ldual(b))).
 * In even dimensions ldual and rdual differ by sign on odd grades. */
#ifndef ALGPY_EXTERIOR2_H
#define ALGPY_EXTERIOR2_H

#include <string.h>

#include "../../vectorspaces/vec4.h"

#define ALGPY_EXTERIOR2_N 4

typedef vec4 exterior2;

static inline exterior2
exterior2_make(const double coeffs[ALGPY_EXTERIOR2_N])
{
    exterior2 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline exterior2
exterior2_grade0(exterior2 a)
{
    const double *x = a.data;
    return exterior2_make((const double[]){x[0], 0, 0, 0});
}

static inline exterior2
exterior2_grade1(exterior2 a)
{
    const double *x = a.data;
    return exterior2_make((const double[]){0, x[1], x[2], 0});
}

static inline exterior2
exterior2_grade2(exterior2 a)
{
    const double *x = a.data;
    return exterior2_make((const double[]){0, 0, 0, x[3]});
}

static inline exterior2
exterior2_from_re(double a)
{
    return exterior2_make((const double[ALGPY_EXTERIOR2_N]){a});
}

static inline double
exterior2_to_re(exterior2 a)
{
    return a.data[0];
}

static inline exterior2
exterior2_zero(void)
{
    return vec4_zero();
}

static inline exterior2
exterior2_bottom(void)
{
    return exterior2_make((const double[ALGPY_EXTERIOR2_N]){[0] = 1});
}

static inline exterior2
exterior2_top(void)
{
    return exterior2_make((const double[ALGPY_EXTERIOR2_N]){[3] = 1});
}

static inline exterior2
exterior2_pos(exterior2 a)
{
    return vec4_pos(a);
}

static inline exterior2
exterior2_neg(exterior2 a)
{
    return vec4_neg(a);
}

static inline exterior2
exterior2_add(exterior2 a, exterior2 b)
{
    return vec4_add(a, b);
}

static inline exterior2
exterior2_sub(exterior2 a, exterior2 b)
{
    return vec4_sub(a, b);
}

static inline exterior2
exterior2_rmul(double a, exterior2 b)
{
    return vec4_scale(b, a);
}

static inline exterior2
exterior2_lmul(exterior2 a, double b)
{
    return vec4_scale(a, b);
}

static inline exterior2
exterior2_ldiv(exterior2 a, double b)
{
    return vec4_div_scalar(a, b);
}

static inline int
exterior2_equal(exterior2 a, exterior2 b)
{
    return vec4_all_eq(a, b);
}

static inline int
exterior2_not_equal(exterior2 a, exterior2 b)
{
    return !exterior2_equal(a, b);
}

static inline int
exterior2_is_real(exterior2 a)
{
    return vec4_is_real(a);
}

static inline exterior2
exterior2_join(exterior2 a, exterior2 b)
{
    const double *x = a.data, *y = b.data;
    return exterior2_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] + x[2]*y[0],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0]});
}

static inline exterior2
exterior2_meet(exterior2 a, exterior2 b)
{
    const double *x = a.data, *y = b.data;
    return exterior2_make((const double[]){
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0],
        x[1]*y[3] + x[3]*y[1],
        x[2]*y[3] + x[3]*y[2],
        x[3]*y[3]});
}

static inline exterior2
exterior2_rdual(exterior2 a)
{
    const double *x = a.data;
    return exterior2_make((const double[]){
        x[3],
        x[2],
        -x[1],
        x[0]});
}

static inline exterior2
exterior2_ldual(exterior2 a)
{
    const double *x = a.data;
    return exterior2_make((const double[]){
        x[3],
        -x[2],
        x[1],
        x[0]});
}

static inline double
exterior2_norm2(exterior2 a)
{
    double result = 0;
    int i;
    for (i = 0; i < ALGPY_EXTERIOR2_N; i++) {
        result += a.data[i]*a.data[i];
    }
    return result;
}

static inline void
exterior2_basis(exterior2 result[4])
{
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
exterior2_gens(exterior2 result[2])
{
    result[0] = exterior2_make((const double[ALGPY_EXTERIOR2_N]){[1] = 1});
    result[1] = exterior2_make((const double[ALGPY_EXTERIOR2_N]){[2] = 1});
}

#endif /* ALGPY_EXTERIOR2_H */
