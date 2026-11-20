/* Exterior algebra on 3 generators as a lattice.
 * The Clifford algebra with every generator square 0; 8 real coefficients.
 * join is the exterior product, with identity bottom = 1.
 * meet is the regressive product, with identity top = e012.
 * rdual(e_a) join e_a = top and e_a join ldual(e_a) = top for basis blades,
 * and meet(a, b) = rdual(join(ldual(a), ldual(b))).
 * In odd dimensions the two duals coincide and are also exposed as dual. */
#ifndef ALGPY_EXTERIOR3_H
#define ALGPY_EXTERIOR3_H

#include <string.h>

#include "../../vectorspaces/vec8.h"

#define ALGPY_EXTERIOR3_N 8

typedef vec8 exterior3;

static inline exterior3
exterior3_make(const double coeffs[ALGPY_EXTERIOR3_N])
{
    exterior3 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline exterior3
exterior3_grade0(exterior3 a)
{
    const double *x = a.data;
    return exterior3_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0});
}

static inline exterior3
exterior3_grade1(exterior3 a)
{
    const double *x = a.data;
    return exterior3_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0});
}

static inline exterior3
exterior3_grade2(exterior3 a)
{
    const double *x = a.data;
    return exterior3_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0});
}

static inline exterior3
exterior3_grade3(exterior3 a)
{
    const double *x = a.data;
    return exterior3_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7]});
}

static inline exterior3
exterior3_from_re(double a)
{
    return exterior3_make((const double[ALGPY_EXTERIOR3_N]){a});
}

static inline double
exterior3_to_re(exterior3 a)
{
    return a.data[0];
}

static inline exterior3
exterior3_zero(void)
{
    return vec8_zero();
}

static inline exterior3
exterior3_bottom(void)
{
    return exterior3_make((const double[ALGPY_EXTERIOR3_N]){[0] = 1});
}

static inline exterior3
exterior3_top(void)
{
    return exterior3_make((const double[ALGPY_EXTERIOR3_N]){[7] = 1});
}

static inline exterior3
exterior3_pos(exterior3 a)
{
    return vec8_pos(a);
}

static inline exterior3
exterior3_neg(exterior3 a)
{
    return vec8_neg(a);
}

static inline exterior3
exterior3_add(exterior3 a, exterior3 b)
{
    return vec8_add(a, b);
}

static inline exterior3
exterior3_sub(exterior3 a, exterior3 b)
{
    return vec8_sub(a, b);
}

static inline exterior3
exterior3_rmul(double a, exterior3 b)
{
    return vec8_scale(b, a);
}

static inline exterior3
exterior3_lmul(exterior3 a, double b)
{
    return vec8_scale(a, b);
}

static inline exterior3
exterior3_ldiv(exterior3 a, double b)
{
    return vec8_div_scalar(a, b);
}

static inline int
exterior3_equal(exterior3 a, exterior3 b)
{
    return vec8_all_eq(a, b);
}

static inline int
exterior3_not_equal(exterior3 a, exterior3 b)
{
    return !exterior3_equal(a, b);
}

static inline int
exterior3_is_real(exterior3 a)
{
    return vec8_is_real(a);
}

static inline exterior3
exterior3_join(exterior3 a, exterior3 b)
{
    const double *x = a.data, *y = b.data;
    return exterior3_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] + x[2]*y[0],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0],
        x[0]*y[4] + x[4]*y[0],
        x[0]*y[5] + x[1]*y[4] - x[4]*y[1] + x[5]*y[0],
        x[0]*y[6] + x[2]*y[4] - x[4]*y[2] + x[6]*y[0],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0]});
}

static inline exterior3
exterior3_meet(exterior3 a, exterior3 b)
{
    const double *x = a.data, *y = b.data;
    return exterior3_make((const double[]){
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0],
        x[1]*y[7] + x[3]*y[5] - x[5]*y[3] + x[7]*y[1],
        x[2]*y[7] + x[3]*y[6] - x[6]*y[3] + x[7]*y[2],
        x[3]*y[7] + x[7]*y[3],
        x[4]*y[7] + x[5]*y[6] - x[6]*y[5] + x[7]*y[4],
        x[5]*y[7] + x[7]*y[5],
        x[6]*y[7] + x[7]*y[6],
        x[7]*y[7]});
}

static inline exterior3
exterior3_rdual(exterior3 a)
{
    const double *x = a.data;
    return exterior3_make((const double[]){
        x[7],
        x[6],
        -x[5],
        x[4],
        x[3],
        -x[2],
        x[1],
        x[0]});
}

static inline exterior3
exterior3_ldual(exterior3 a)
{
    const double *x = a.data;
    return exterior3_make((const double[]){
        x[7],
        x[6],
        -x[5],
        x[4],
        x[3],
        -x[2],
        x[1],
        x[0]});
}

static inline exterior3
exterior3_dual(exterior3 a)
{
    return exterior3_rdual(a);
}

static inline double
exterior3_norm2(exterior3 a)
{
    double result = 0;
    int i;
    for (i = 0; i < ALGPY_EXTERIOR3_N; i++) {
        result += a.data[i]*a.data[i];
    }
    return result;
}

static inline void
exterior3_basis(exterior3 result[8])
{
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
exterior3_gens(exterior3 result[3])
{
    result[0] = exterior3_make((const double[ALGPY_EXTERIOR3_N]){[1] = 1});
    result[1] = exterior3_make((const double[ALGPY_EXTERIOR3_N]){[2] = 1});
    result[2] = exterior3_make((const double[ALGPY_EXTERIOR3_N]){[4] = 1});
}

#endif /* ALGPY_EXTERIOR3_H */
