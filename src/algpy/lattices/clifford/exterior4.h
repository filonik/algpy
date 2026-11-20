/* Exterior algebra on 4 generators as a lattice.
 * The Clifford algebra with every generator square 0; 16 real coefficients.
 * join is the exterior product, with identity bottom = 1.
 * meet is the regressive product, with identity top = e0123.
 * rdual(e_a) join e_a = top and e_a join ldual(e_a) = top for basis blades,
 * and meet(a, b) = rdual(join(ldual(a), ldual(b))).
 * In even dimensions ldual and rdual differ by sign on odd grades. */
#ifndef ALGPY_EXTERIOR4_H
#define ALGPY_EXTERIOR4_H

#include <string.h>

#include "../../vectorspaces/vec16.h"

#define ALGPY_EXTERIOR4_N 16

typedef vec16 exterior4;

static inline exterior4
exterior4_make(const double coeffs[ALGPY_EXTERIOR4_N])
{
    exterior4 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline exterior4
exterior4_grade0(exterior4 a)
{
    const double *x = a.data;
    return exterior4_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline exterior4
exterior4_grade1(exterior4 a)
{
    const double *x = a.data;
    return exterior4_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0});
}

static inline exterior4
exterior4_grade2(exterior4 a)
{
    const double *x = a.data;
    return exterior4_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0});
}

static inline exterior4
exterior4_grade3(exterior4 a)
{
    const double *x = a.data;
    return exterior4_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0});
}

static inline exterior4
exterior4_grade4(exterior4 a)
{
    const double *x = a.data;
    return exterior4_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15]});
}

static inline exterior4
exterior4_from_re(double a)
{
    return exterior4_make((const double[ALGPY_EXTERIOR4_N]){a});
}

static inline double
exterior4_to_re(exterior4 a)
{
    return a.data[0];
}

static inline exterior4
exterior4_zero(void)
{
    return vec16_zero();
}

static inline exterior4
exterior4_bottom(void)
{
    return exterior4_make((const double[ALGPY_EXTERIOR4_N]){[0] = 1});
}

static inline exterior4
exterior4_top(void)
{
    return exterior4_make((const double[ALGPY_EXTERIOR4_N]){[15] = 1});
}

static inline exterior4
exterior4_pos(exterior4 a)
{
    return vec16_pos(a);
}

static inline exterior4
exterior4_neg(exterior4 a)
{
    return vec16_neg(a);
}

static inline exterior4
exterior4_add(exterior4 a, exterior4 b)
{
    return vec16_add(a, b);
}

static inline exterior4
exterior4_sub(exterior4 a, exterior4 b)
{
    return vec16_sub(a, b);
}

static inline exterior4
exterior4_rmul(double a, exterior4 b)
{
    return vec16_scale(b, a);
}

static inline exterior4
exterior4_lmul(exterior4 a, double b)
{
    return vec16_scale(a, b);
}

static inline exterior4
exterior4_ldiv(exterior4 a, double b)
{
    return vec16_div_scalar(a, b);
}

static inline int
exterior4_equal(exterior4 a, exterior4 b)
{
    return vec16_all_eq(a, b);
}

static inline int
exterior4_not_equal(exterior4 a, exterior4 b)
{
    return !exterior4_equal(a, b);
}

static inline int
exterior4_is_real(exterior4 a)
{
    return vec16_is_real(a);
}

static inline exterior4
exterior4_join(exterior4 a, exterior4 b)
{
    const double *x = a.data, *y = b.data;
    return exterior4_make((const double[]){
        x[0]*y[0],
        x[0]*y[1] + x[1]*y[0],
        x[0]*y[2] + x[2]*y[0],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0],
        x[0]*y[4] + x[4]*y[0],
        x[0]*y[5] + x[1]*y[4] - x[4]*y[1] + x[5]*y[0],
        x[0]*y[6] + x[2]*y[4] - x[4]*y[2] + x[6]*y[0],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0],
        x[0]*y[8] + x[8]*y[0],
        x[0]*y[9] + x[1]*y[8] - x[8]*y[1] + x[9]*y[0],
        x[0]*y[10] + x[2]*y[8] - x[8]*y[2] + x[10]*y[0],
        x[0]*y[11] + x[1]*y[10] - x[2]*y[9] + x[3]*y[8] + x[8]*y[3] - x[9]*y[2] + x[10]*y[1] + x[11]*y[0],
        x[0]*y[12] + x[4]*y[8] - x[8]*y[4] + x[12]*y[0],
        x[0]*y[13] + x[1]*y[12] - x[4]*y[9] + x[5]*y[8] + x[8]*y[5] - x[9]*y[4] + x[12]*y[1] + x[13]*y[0],
        x[0]*y[14] + x[2]*y[12] - x[4]*y[10] + x[6]*y[8] + x[8]*y[6] - x[10]*y[4] + x[12]*y[2] + x[14]*y[0],
        x[0]*y[15] + x[1]*y[14] - x[2]*y[13] + x[3]*y[12] + x[4]*y[11] - x[5]*y[10] + x[6]*y[9] + x[7]*y[8] - x[8]*y[7] + x[9]*y[6] - x[10]*y[5] - x[11]*y[4] + x[12]*y[3] + x[13]*y[2] - x[14]*y[1] + x[15]*y[0]});
}

static inline exterior4
exterior4_meet(exterior4 a, exterior4 b)
{
    const double *x = a.data, *y = b.data;
    return exterior4_make((const double[]){
        x[0]*y[15] + x[1]*y[14] - x[2]*y[13] + x[3]*y[12] + x[4]*y[11] - x[5]*y[10] + x[6]*y[9] + x[7]*y[8] - x[8]*y[7] + x[9]*y[6] - x[10]*y[5] - x[11]*y[4] + x[12]*y[3] + x[13]*y[2] - x[14]*y[1] + x[15]*y[0],
        x[1]*y[15] + x[3]*y[13] - x[5]*y[11] + x[7]*y[9] + x[9]*y[7] - x[11]*y[5] + x[13]*y[3] + x[15]*y[1],
        x[2]*y[15] + x[3]*y[14] - x[6]*y[11] + x[7]*y[10] + x[10]*y[7] - x[11]*y[6] + x[14]*y[3] + x[15]*y[2],
        x[3]*y[15] + x[7]*y[11] - x[11]*y[7] + x[15]*y[3],
        x[4]*y[15] + x[5]*y[14] - x[6]*y[13] + x[7]*y[12] + x[12]*y[7] - x[13]*y[6] + x[14]*y[5] + x[15]*y[4],
        x[5]*y[15] + x[7]*y[13] - x[13]*y[7] + x[15]*y[5],
        x[6]*y[15] + x[7]*y[14] - x[14]*y[7] + x[15]*y[6],
        x[7]*y[15] + x[15]*y[7],
        x[8]*y[15] + x[9]*y[14] - x[10]*y[13] + x[11]*y[12] + x[12]*y[11] - x[13]*y[10] + x[14]*y[9] + x[15]*y[8],
        x[9]*y[15] + x[11]*y[13] - x[13]*y[11] + x[15]*y[9],
        x[10]*y[15] + x[11]*y[14] - x[14]*y[11] + x[15]*y[10],
        x[11]*y[15] + x[15]*y[11],
        x[12]*y[15] + x[13]*y[14] - x[14]*y[13] + x[15]*y[12],
        x[13]*y[15] + x[15]*y[13],
        x[14]*y[15] + x[15]*y[14],
        x[15]*y[15]});
}

static inline exterior4
exterior4_rdual(exterior4 a)
{
    const double *x = a.data;
    return exterior4_make((const double[]){
        x[15],
        x[14],
        -x[13],
        x[12],
        x[11],
        -x[10],
        x[9],
        x[8],
        -x[7],
        x[6],
        -x[5],
        -x[4],
        x[3],
        x[2],
        -x[1],
        x[0]});
}

static inline exterior4
exterior4_ldual(exterior4 a)
{
    const double *x = a.data;
    return exterior4_make((const double[]){
        x[15],
        -x[14],
        x[13],
        x[12],
        -x[11],
        -x[10],
        x[9],
        -x[8],
        x[7],
        x[6],
        -x[5],
        x[4],
        x[3],
        -x[2],
        x[1],
        x[0]});
}

static inline double
exterior4_norm2(exterior4 a)
{
    double result = 0;
    int i;
    for (i = 0; i < ALGPY_EXTERIOR4_N; i++) {
        result += a.data[i]*a.data[i];
    }
    return result;
}

static inline void
exterior4_basis(exterior4 result[16])
{
    int i, j;
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
exterior4_gens(exterior4 result[4])
{
    result[0] = exterior4_make((const double[ALGPY_EXTERIOR4_N]){[1] = 1});
    result[1] = exterior4_make((const double[ALGPY_EXTERIOR4_N]){[2] = 1});
    result[2] = exterior4_make((const double[ALGPY_EXTERIOR4_N]){[4] = 1});
    result[3] = exterior4_make((const double[ALGPY_EXTERIOR4_N]){[8] = 1});
}

#endif /* ALGPY_EXTERIOR4_H */
