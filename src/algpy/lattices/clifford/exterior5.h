/* Exterior algebra on 5 generators as a lattice.
 * The Clifford algebra with every generator square 0; 32 real coefficients.
 * join is the exterior product, with identity bottom = 1.
 * meet is the regressive product, with identity top = e01234.
 * rdual(e_a) join e_a = top and e_a join ldual(e_a) = top for basis blades,
 * and meet(a, b) = rdual(join(ldual(a), ldual(b))).
 * In odd dimensions the two duals coincide and are also exposed as dual. */
#ifndef ALGPY_EXTERIOR5_H
#define ALGPY_EXTERIOR5_H

#include <string.h>

#include "../../vectorspaces/vec32.h"

#define ALGPY_EXTERIOR5_N 32

typedef vec32 exterior5;

static inline exterior5
exterior5_make(const double coeffs[ALGPY_EXTERIOR5_N])
{
    exterior5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline exterior5
exterior5_grade0(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline exterior5
exterior5_grade1(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline exterior5
exterior5_grade2(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline exterior5
exterior5_grade3(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline exterior5
exterior5_grade4(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline exterior5
exterior5_grade5(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

static inline exterior5
exterior5_from_re(double a)
{
    return exterior5_make((const double[ALGPY_EXTERIOR5_N]){a});
}

static inline double
exterior5_to_re(exterior5 a)
{
    return a.data[0];
}

static inline exterior5
exterior5_zero(void)
{
    return vec32_zero();
}

static inline exterior5
exterior5_bottom(void)
{
    return exterior5_make((const double[ALGPY_EXTERIOR5_N]){[0] = 1});
}

static inline exterior5
exterior5_top(void)
{
    return exterior5_make((const double[ALGPY_EXTERIOR5_N]){[31] = 1});
}

static inline exterior5
exterior5_pos(exterior5 a)
{
    return vec32_pos(a);
}

static inline exterior5
exterior5_neg(exterior5 a)
{
    return vec32_neg(a);
}

static inline exterior5
exterior5_add(exterior5 a, exterior5 b)
{
    return vec32_add(a, b);
}

static inline exterior5
exterior5_sub(exterior5 a, exterior5 b)
{
    return vec32_sub(a, b);
}

static inline exterior5
exterior5_rmul(double a, exterior5 b)
{
    return vec32_scale(b, a);
}

static inline exterior5
exterior5_lmul(exterior5 a, double b)
{
    return vec32_scale(a, b);
}

static inline exterior5
exterior5_ldiv(exterior5 a, double b)
{
    return vec32_div_scalar(a, b);
}

static inline int
exterior5_equal(exterior5 a, exterior5 b)
{
    return vec32_all_eq(a, b);
}

static inline int
exterior5_not_equal(exterior5 a, exterior5 b)
{
    return !exterior5_equal(a, b);
}

static inline int
exterior5_is_real(exterior5 a)
{
    return vec32_is_real(a);
}

static inline exterior5
exterior5_join(exterior5 a, exterior5 b)
{
    const double *x = a.data, *y = b.data;
    return exterior5_make((const double[]){
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
        x[0]*y[15] + x[1]*y[14] - x[2]*y[13] + x[3]*y[12] + x[4]*y[11] - x[5]*y[10] + x[6]*y[9] + x[7]*y[8] - x[8]*y[7] + x[9]*y[6] - x[10]*y[5] - x[11]*y[4] + x[12]*y[3] + x[13]*y[2] - x[14]*y[1] + x[15]*y[0],
        x[0]*y[16] + x[16]*y[0],
        x[0]*y[17] + x[1]*y[16] - x[16]*y[1] + x[17]*y[0],
        x[0]*y[18] + x[2]*y[16] - x[16]*y[2] + x[18]*y[0],
        x[0]*y[19] + x[1]*y[18] - x[2]*y[17] + x[3]*y[16] + x[16]*y[3] - x[17]*y[2] + x[18]*y[1] + x[19]*y[0],
        x[0]*y[20] + x[4]*y[16] - x[16]*y[4] + x[20]*y[0],
        x[0]*y[21] + x[1]*y[20] - x[4]*y[17] + x[5]*y[16] + x[16]*y[5] - x[17]*y[4] + x[20]*y[1] + x[21]*y[0],
        x[0]*y[22] + x[2]*y[20] - x[4]*y[18] + x[6]*y[16] + x[16]*y[6] - x[18]*y[4] + x[20]*y[2] + x[22]*y[0],
        x[0]*y[23] + x[1]*y[22] - x[2]*y[21] + x[3]*y[20] + x[4]*y[19] - x[5]*y[18] + x[6]*y[17] + x[7]*y[16] - x[16]*y[7] + x[17]*y[6] - x[18]*y[5] - x[19]*y[4] + x[20]*y[3] + x[21]*y[2] - x[22]*y[1] + x[23]*y[0],
        x[0]*y[24] + x[8]*y[16] - x[16]*y[8] + x[24]*y[0],
        x[0]*y[25] + x[1]*y[24] - x[8]*y[17] + x[9]*y[16] + x[16]*y[9] - x[17]*y[8] + x[24]*y[1] + x[25]*y[0],
        x[0]*y[26] + x[2]*y[24] - x[8]*y[18] + x[10]*y[16] + x[16]*y[10] - x[18]*y[8] + x[24]*y[2] + x[26]*y[0],
        x[0]*y[27] + x[1]*y[26] - x[2]*y[25] + x[3]*y[24] + x[8]*y[19] - x[9]*y[18] + x[10]*y[17] + x[11]*y[16] - x[16]*y[11] + x[17]*y[10] - x[18]*y[9] - x[19]*y[8] + x[24]*y[3] + x[25]*y[2] - x[26]*y[1] + x[27]*y[0],
        x[0]*y[28] + x[4]*y[24] - x[8]*y[20] + x[12]*y[16] + x[16]*y[12] - x[20]*y[8] + x[24]*y[4] + x[28]*y[0],
        x[0]*y[29] + x[1]*y[28] - x[4]*y[25] + x[5]*y[24] + x[8]*y[21] - x[9]*y[20] + x[12]*y[17] + x[13]*y[16] - x[16]*y[13] + x[17]*y[12] - x[20]*y[9] - x[21]*y[8] + x[24]*y[5] + x[25]*y[4] - x[28]*y[1] + x[29]*y[0],
        x[0]*y[30] + x[2]*y[28] - x[4]*y[26] + x[6]*y[24] + x[8]*y[22] - x[10]*y[20] + x[12]*y[18] + x[14]*y[16] - x[16]*y[14] + x[18]*y[12] - x[20]*y[10] - x[22]*y[8] + x[24]*y[6] + x[26]*y[4] - x[28]*y[2] + x[30]*y[0],
        x[0]*y[31] + x[1]*y[30] - x[2]*y[29] + x[3]*y[28] + x[4]*y[27] - x[5]*y[26] + x[6]*y[25] + x[7]*y[24] - x[8]*y[23] + x[9]*y[22] - x[10]*y[21] - x[11]*y[20] + x[12]*y[19] + x[13]*y[18] - x[14]*y[17] + x[15]*y[16] + x[16]*y[15] - x[17]*y[14] + x[18]*y[13] + x[19]*y[12] - x[20]*y[11] - x[21]*y[10] + x[22]*y[9] - x[23]*y[8] + x[24]*y[7] + x[25]*y[6] - x[26]*y[5] + x[27]*y[4] + x[28]*y[3] - x[29]*y[2] + x[30]*y[1] + x[31]*y[0]});
}

static inline exterior5
exterior5_meet(exterior5 a, exterior5 b)
{
    const double *x = a.data, *y = b.data;
    return exterior5_make((const double[]){
        x[0]*y[31] + x[1]*y[30] - x[2]*y[29] + x[3]*y[28] + x[4]*y[27] - x[5]*y[26] + x[6]*y[25] + x[7]*y[24] - x[8]*y[23] + x[9]*y[22] - x[10]*y[21] - x[11]*y[20] + x[12]*y[19] + x[13]*y[18] - x[14]*y[17] + x[15]*y[16] + x[16]*y[15] - x[17]*y[14] + x[18]*y[13] + x[19]*y[12] - x[20]*y[11] - x[21]*y[10] + x[22]*y[9] - x[23]*y[8] + x[24]*y[7] + x[25]*y[6] - x[26]*y[5] + x[27]*y[4] + x[28]*y[3] - x[29]*y[2] + x[30]*y[1] + x[31]*y[0],
        x[1]*y[31] + x[3]*y[29] - x[5]*y[27] + x[7]*y[25] + x[9]*y[23] - x[11]*y[21] + x[13]*y[19] + x[15]*y[17] - x[17]*y[15] + x[19]*y[13] - x[21]*y[11] - x[23]*y[9] + x[25]*y[7] + x[27]*y[5] - x[29]*y[3] + x[31]*y[1],
        x[2]*y[31] + x[3]*y[30] - x[6]*y[27] + x[7]*y[26] + x[10]*y[23] - x[11]*y[22] + x[14]*y[19] + x[15]*y[18] - x[18]*y[15] + x[19]*y[14] - x[22]*y[11] - x[23]*y[10] + x[26]*y[7] + x[27]*y[6] - x[30]*y[3] + x[31]*y[2],
        x[3]*y[31] + x[7]*y[27] - x[11]*y[23] + x[15]*y[19] + x[19]*y[15] - x[23]*y[11] + x[27]*y[7] + x[31]*y[3],
        x[4]*y[31] + x[5]*y[30] - x[6]*y[29] + x[7]*y[28] + x[12]*y[23] - x[13]*y[22] + x[14]*y[21] + x[15]*y[20] - x[20]*y[15] + x[21]*y[14] - x[22]*y[13] - x[23]*y[12] + x[28]*y[7] + x[29]*y[6] - x[30]*y[5] + x[31]*y[4],
        x[5]*y[31] + x[7]*y[29] - x[13]*y[23] + x[15]*y[21] + x[21]*y[15] - x[23]*y[13] + x[29]*y[7] + x[31]*y[5],
        x[6]*y[31] + x[7]*y[30] - x[14]*y[23] + x[15]*y[22] + x[22]*y[15] - x[23]*y[14] + x[30]*y[7] + x[31]*y[6],
        x[7]*y[31] + x[15]*y[23] - x[23]*y[15] + x[31]*y[7],
        x[8]*y[31] + x[9]*y[30] - x[10]*y[29] + x[11]*y[28] + x[12]*y[27] - x[13]*y[26] + x[14]*y[25] + x[15]*y[24] - x[24]*y[15] + x[25]*y[14] - x[26]*y[13] - x[27]*y[12] + x[28]*y[11] + x[29]*y[10] - x[30]*y[9] + x[31]*y[8],
        x[9]*y[31] + x[11]*y[29] - x[13]*y[27] + x[15]*y[25] + x[25]*y[15] - x[27]*y[13] + x[29]*y[11] + x[31]*y[9],
        x[10]*y[31] + x[11]*y[30] - x[14]*y[27] + x[15]*y[26] + x[26]*y[15] - x[27]*y[14] + x[30]*y[11] + x[31]*y[10],
        x[11]*y[31] + x[15]*y[27] - x[27]*y[15] + x[31]*y[11],
        x[12]*y[31] + x[13]*y[30] - x[14]*y[29] + x[15]*y[28] + x[28]*y[15] - x[29]*y[14] + x[30]*y[13] + x[31]*y[12],
        x[13]*y[31] + x[15]*y[29] - x[29]*y[15] + x[31]*y[13],
        x[14]*y[31] + x[15]*y[30] - x[30]*y[15] + x[31]*y[14],
        x[15]*y[31] + x[31]*y[15],
        x[16]*y[31] + x[17]*y[30] - x[18]*y[29] + x[19]*y[28] + x[20]*y[27] - x[21]*y[26] + x[22]*y[25] + x[23]*y[24] - x[24]*y[23] + x[25]*y[22] - x[26]*y[21] - x[27]*y[20] + x[28]*y[19] + x[29]*y[18] - x[30]*y[17] + x[31]*y[16],
        x[17]*y[31] + x[19]*y[29] - x[21]*y[27] + x[23]*y[25] + x[25]*y[23] - x[27]*y[21] + x[29]*y[19] + x[31]*y[17],
        x[18]*y[31] + x[19]*y[30] - x[22]*y[27] + x[23]*y[26] + x[26]*y[23] - x[27]*y[22] + x[30]*y[19] + x[31]*y[18],
        x[19]*y[31] + x[23]*y[27] - x[27]*y[23] + x[31]*y[19],
        x[20]*y[31] + x[21]*y[30] - x[22]*y[29] + x[23]*y[28] + x[28]*y[23] - x[29]*y[22] + x[30]*y[21] + x[31]*y[20],
        x[21]*y[31] + x[23]*y[29] - x[29]*y[23] + x[31]*y[21],
        x[22]*y[31] + x[23]*y[30] - x[30]*y[23] + x[31]*y[22],
        x[23]*y[31] + x[31]*y[23],
        x[24]*y[31] + x[25]*y[30] - x[26]*y[29] + x[27]*y[28] + x[28]*y[27] - x[29]*y[26] + x[30]*y[25] + x[31]*y[24],
        x[25]*y[31] + x[27]*y[29] - x[29]*y[27] + x[31]*y[25],
        x[26]*y[31] + x[27]*y[30] - x[30]*y[27] + x[31]*y[26],
        x[27]*y[31] + x[31]*y[27],
        x[28]*y[31] + x[29]*y[30] - x[30]*y[29] + x[31]*y[28],
        x[29]*y[31] + x[31]*y[29],
        x[30]*y[31] + x[31]*y[30],
        x[31]*y[31]});
}

static inline exterior5
exterior5_rdual(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){
        x[31],
        x[30],
        -x[29],
        x[28],
        x[27],
        -x[26],
        x[25],
        x[24],
        -x[23],
        x[22],
        -x[21],
        -x[20],
        x[19],
        x[18],
        -x[17],
        x[16],
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

static inline exterior5
exterior5_ldual(exterior5 a)
{
    const double *x = a.data;
    return exterior5_make((const double[]){
        x[31],
        x[30],
        -x[29],
        x[28],
        x[27],
        -x[26],
        x[25],
        x[24],
        -x[23],
        x[22],
        -x[21],
        -x[20],
        x[19],
        x[18],
        -x[17],
        x[16],
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

static inline exterior5
exterior5_dual(exterior5 a)
{
    return exterior5_rdual(a);
}

static inline double
exterior5_norm2(exterior5 a)
{
    double result = 0;
    int i;
    for (i = 0; i < ALGPY_EXTERIOR5_N; i++) {
        result += a.data[i]*a.data[i];
    }
    return result;
}

static inline void
exterior5_basis(exterior5 result[32])
{
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
exterior5_gens(exterior5 result[5])
{
    result[0] = exterior5_make((const double[ALGPY_EXTERIOR5_N]){[1] = 1});
    result[1] = exterior5_make((const double[ALGPY_EXTERIOR5_N]){[2] = 1});
    result[2] = exterior5_make((const double[ALGPY_EXTERIOR5_N]){[4] = 1});
    result[3] = exterior5_make((const double[ALGPY_EXTERIOR5_N]){[8] = 1});
    result[4] = exterior5_make((const double[ALGPY_EXTERIOR5_N]){[16] = 1});
}

#endif /* ALGPY_EXTERIOR5_H */
