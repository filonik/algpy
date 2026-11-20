/* Conformal geometric algebra in 3D.
 * Generator squares: -1, 1, 1, 1, 1.
 * 5 generators; 32 real coefficients. */
#ifndef ALGPY_CGA3_H
#define ALGPY_CGA3_H

#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../polynomials/hermite_complex.h"
#include "../../vectorspaces/vec32.h"
#include "../../lattices/clifford/exterior5.h"
#include "../multicomplex/complex1.h"
#include "cga2.h"

#define ALGPY_CGA3_N 32

typedef vec32 cga3;

static inline cga3
cga3_make(const double coeffs[ALGPY_CGA3_N])
{
    cga3 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline cga3
cga3_grade0(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline cga3
cga3_grade1(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline cga3
cga3_grade2(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline cga3
cga3_grade3(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline cga3
cga3_grade4(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline cga3
cga3_grade5(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

static inline cga3
cga3_from_re(double a)
{
    cga3 result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
cga3_to_re(cga3 a)
{
    return a.data[0];
}

static inline cga3
cga3_zero(void)
{
    return vec32_zero();
}

static inline cga3
cga3_one(void)
{
    return cga3_from_re(1);
}

static inline cga3
cga3_pos(cga3 a)
{
    return vec32_pos(a);
}

static inline cga3
cga3_neg(cga3 a)
{
    return vec32_neg(a);
}

static inline cga3
cga3_add(cga3 a, cga3 b)
{
    return vec32_add(a, b);
}

static inline cga3
cga3_sub(cga3 a, cga3 b)
{
    return vec32_sub(a, b);
}

static inline int
cga3_equal(cga3 a, cga3 b)
{
    return vec32_all_eq(a, b);
}

static inline int
cga3_not_equal(cga3 a, cga3 b)
{
    return !cga3_equal(a, b);
}

static inline int
cga3_is_real(cga3 a)
{
    return vec32_is_real(a);
}

static inline cga3
cga3_re(cga3 a)
{
    return cga3_from_re(a.data[0]);
}

static inline cga3
cga3_im(cga3 a)
{
    return cga3_sub(a, cga3_re(a));
}

static inline cga3
cga3_mul(cga3 a, cga3 b)
{
    const double *x = a.data, *y = b.data;
    return cga3_make((const double[]){
        x[0]*y[0] - x[1]*y[1] + x[2]*y[2] + x[3]*y[3] + x[4]*y[4] + x[5]*y[5] - x[6]*y[6] + x[7]*y[7] + x[8]*y[8] + x[9]*y[9] - x[10]*y[10] + x[11]*y[11] - x[12]*y[12] + x[13]*y[13] - x[14]*y[14] - x[15]*y[15] + x[16]*y[16] + x[17]*y[17] - x[18]*y[18] + x[19]*y[19] - x[20]*y[20] + x[21]*y[21] - x[22]*y[22] - x[23]*y[23] - x[24]*y[24] + x[25]*y[25] - x[26]*y[26] - x[27]*y[27] - x[28]*y[28] - x[29]*y[29] + x[30]*y[30] - x[31]*y[31],
        x[0]*y[1] + x[1]*y[0] - x[2]*y[3] + x[3]*y[2] - x[4]*y[5] + x[5]*y[4] - x[6]*y[7] - x[7]*y[6] - x[8]*y[9] + x[9]*y[8] - x[10]*y[11] - x[11]*y[10] - x[12]*y[13] - x[13]*y[12] + x[14]*y[15] - x[15]*y[14] - x[16]*y[17] + x[17]*y[16] - x[18]*y[19] - x[19]*y[18] - x[20]*y[21] - x[21]*y[20] + x[22]*y[23] - x[23]*y[22] - x[24]*y[25] - x[25]*y[24] + x[26]*y[27] - x[27]*y[26] + x[28]*y[29] - x[29]*y[28] + x[30]*y[31] + x[31]*y[30],
        x[0]*y[2] - x[1]*y[3] + x[2]*y[0] + x[3]*y[1] - x[4]*y[6] - x[5]*y[7] + x[6]*y[4] - x[7]*y[5] - x[8]*y[10] - x[9]*y[11] + x[10]*y[8] - x[11]*y[9] - x[12]*y[14] + x[13]*y[15] - x[14]*y[12] - x[15]*y[13] - x[16]*y[18] - x[17]*y[19] + x[18]*y[16] - x[19]*y[17] - x[20]*y[22] + x[21]*y[23] - x[22]*y[20] - x[23]*y[21] - x[24]*y[26] + x[25]*y[27] - x[26]*y[24] - x[27]*y[25] + x[28]*y[30] + x[29]*y[31] - x[30]*y[28] + x[31]*y[29],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0] + x[4]*y[7] - x[5]*y[6] + x[6]*y[5] + x[7]*y[4] + x[8]*y[11] - x[9]*y[10] + x[10]*y[9] + x[11]*y[8] - x[12]*y[15] - x[13]*y[14] + x[14]*y[13] - x[15]*y[12] + x[16]*y[19] - x[17]*y[18] + x[18]*y[17] + x[19]*y[16] - x[20]*y[23] - x[21]*y[22] + x[22]*y[21] - x[23]*y[20] - x[24]*y[27] - x[25]*y[26] + x[26]*y[25] - x[27]*y[24] - x[28]*y[31] + x[29]*y[30] - x[30]*y[29] - x[31]*y[28],
        x[0]*y[4] - x[1]*y[5] + x[2]*y[6] + x[3]*y[7] + x[4]*y[0] + x[5]*y[1] - x[6]*y[2] + x[7]*y[3] - x[8]*y[12] - x[9]*y[13] + x[10]*y[14] - x[11]*y[15] + x[12]*y[8] - x[13]*y[9] + x[14]*y[10] + x[15]*y[11] - x[16]*y[20] - x[17]*y[21] + x[18]*y[22] - x[19]*y[23] + x[20]*y[16] - x[21]*y[17] + x[22]*y[18] + x[23]*y[19] - x[24]*y[28] + x[25]*y[29] - x[26]*y[30] - x[27]*y[31] - x[28]*y[24] - x[29]*y[25] + x[30]*y[26] - x[31]*y[27],
        x[0]*y[5] + x[1]*y[4] - x[2]*y[7] + x[3]*y[6] - x[4]*y[1] + x[5]*y[0] - x[6]*y[3] - x[7]*y[2] + x[8]*y[13] - x[9]*y[12] + x[10]*y[15] + x[11]*y[14] + x[12]*y[9] + x[13]*y[8] - x[14]*y[11] + x[15]*y[10] + x[16]*y[21] - x[17]*y[20] + x[18]*y[23] + x[19]*y[22] + x[20]*y[17] + x[21]*y[16] - x[22]*y[19] + x[23]*y[18] - x[24]*y[29] - x[25]*y[28] + x[26]*y[31] - x[27]*y[30] + x[28]*y[25] - x[29]*y[24] + x[30]*y[27] + x[31]*y[26],
        x[0]*y[6] - x[1]*y[7] + x[2]*y[4] + x[3]*y[5] - x[4]*y[2] - x[5]*y[3] + x[6]*y[0] - x[7]*y[1] + x[8]*y[14] + x[9]*y[15] - x[10]*y[12] + x[11]*y[13] + x[12]*y[10] - x[13]*y[11] + x[14]*y[8] + x[15]*y[9] + x[16]*y[22] + x[17]*y[23] - x[18]*y[20] + x[19]*y[21] + x[20]*y[18] - x[21]*y[19] + x[22]*y[16] + x[23]*y[17] - x[24]*y[30] + x[25]*y[31] - x[26]*y[28] - x[27]*y[29] + x[28]*y[26] + x[29]*y[27] - x[30]*y[24] + x[31]*y[25],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0] - x[8]*y[15] + x[9]*y[14] - x[10]*y[13] - x[11]*y[12] + x[12]*y[11] + x[13]*y[10] - x[14]*y[9] + x[15]*y[8] - x[16]*y[23] + x[17]*y[22] - x[18]*y[21] - x[19]*y[20] + x[20]*y[19] + x[21]*y[18] - x[22]*y[17] + x[23]*y[16] - x[24]*y[31] - x[25]*y[30] + x[26]*y[29] - x[27]*y[28] - x[28]*y[27] + x[29]*y[26] - x[30]*y[25] - x[31]*y[24],
        x[0]*y[8] - x[1]*y[9] + x[2]*y[10] + x[3]*y[11] + x[4]*y[12] + x[5]*y[13] - x[6]*y[14] + x[7]*y[15] + x[8]*y[0] + x[9]*y[1] - x[10]*y[2] + x[11]*y[3] - x[12]*y[4] + x[13]*y[5] - x[14]*y[6] - x[15]*y[7] - x[16]*y[24] - x[17]*y[25] + x[18]*y[26] - x[19]*y[27] + x[20]*y[28] - x[21]*y[29] + x[22]*y[30] + x[23]*y[31] + x[24]*y[16] - x[25]*y[17] + x[26]*y[18] + x[27]*y[19] + x[28]*y[20] + x[29]*y[21] - x[30]*y[22] + x[31]*y[23],
        x[0]*y[9] + x[1]*y[8] - x[2]*y[11] + x[3]*y[10] - x[4]*y[13] + x[5]*y[12] - x[6]*y[15] - x[7]*y[14] - x[8]*y[1] + x[9]*y[0] - x[10]*y[3] - x[11]*y[2] - x[12]*y[5] - x[13]*y[4] + x[14]*y[7] - x[15]*y[6] + x[16]*y[25] - x[17]*y[24] + x[18]*y[27] + x[19]*y[26] + x[20]*y[29] + x[21]*y[28] - x[22]*y[31] + x[23]*y[30] + x[24]*y[17] + x[25]*y[16] - x[26]*y[19] + x[27]*y[18] - x[28]*y[21] + x[29]*y[20] - x[30]*y[23] - x[31]*y[22],
        x[0]*y[10] - x[1]*y[11] + x[2]*y[8] + x[3]*y[9] - x[4]*y[14] - x[5]*y[15] + x[6]*y[12] - x[7]*y[13] - x[8]*y[2] - x[9]*y[3] + x[10]*y[0] - x[11]*y[1] - x[12]*y[6] + x[13]*y[7] - x[14]*y[4] - x[15]*y[5] + x[16]*y[26] + x[17]*y[27] - x[18]*y[24] + x[19]*y[25] + x[20]*y[30] - x[21]*y[31] + x[22]*y[28] + x[23]*y[29] + x[24]*y[18] - x[25]*y[19] + x[26]*y[16] + x[27]*y[17] - x[28]*y[22] - x[29]*y[23] + x[30]*y[20] - x[31]*y[21],
        x[0]*y[11] + x[1]*y[10] - x[2]*y[9] + x[3]*y[8] + x[4]*y[15] - x[5]*y[14] + x[6]*y[13] + x[7]*y[12] + x[8]*y[3] - x[9]*y[2] + x[10]*y[1] + x[11]*y[0] - x[12]*y[7] - x[13]*y[6] + x[14]*y[5] - x[15]*y[4] - x[16]*y[27] + x[17]*y[26] - x[18]*y[25] - x[19]*y[24] + x[20]*y[31] + x[21]*y[30] - x[22]*y[29] + x[23]*y[28] + x[24]*y[19] + x[25]*y[18] - x[26]*y[17] + x[27]*y[16] + x[28]*y[23] - x[29]*y[22] + x[30]*y[21] + x[31]*y[20],
        x[0]*y[12] - x[1]*y[13] + x[2]*y[14] + x[3]*y[15] + x[4]*y[8] + x[5]*y[9] - x[6]*y[10] + x[7]*y[11] - x[8]*y[4] - x[9]*y[5] + x[10]*y[6] - x[11]*y[7] + x[12]*y[0] - x[13]*y[1] + x[14]*y[2] + x[15]*y[3] + x[16]*y[28] + x[17]*y[29] - x[18]*y[30] + x[19]*y[31] - x[20]*y[24] + x[21]*y[25] - x[22]*y[26] - x[23]*y[27] + x[24]*y[20] - x[25]*y[21] + x[26]*y[22] + x[27]*y[23] + x[28]*y[16] + x[29]*y[17] - x[30]*y[18] + x[31]*y[19],
        x[0]*y[13] + x[1]*y[12] - x[2]*y[15] + x[3]*y[14] - x[4]*y[9] + x[5]*y[8] - x[6]*y[11] - x[7]*y[10] + x[8]*y[5] - x[9]*y[4] + x[10]*y[7] + x[11]*y[6] + x[12]*y[1] + x[13]*y[0] - x[14]*y[3] + x[15]*y[2] - x[16]*y[29] + x[17]*y[28] - x[18]*y[31] - x[19]*y[30] - x[20]*y[25] - x[21]*y[24] + x[22]*y[27] - x[23]*y[26] + x[24]*y[21] + x[25]*y[20] - x[26]*y[23] + x[27]*y[22] - x[28]*y[17] + x[29]*y[16] - x[30]*y[19] - x[31]*y[18],
        x[0]*y[14] - x[1]*y[15] + x[2]*y[12] + x[3]*y[13] - x[4]*y[10] - x[5]*y[11] + x[6]*y[8] - x[7]*y[9] + x[8]*y[6] + x[9]*y[7] - x[10]*y[4] + x[11]*y[5] + x[12]*y[2] - x[13]*y[3] + x[14]*y[0] + x[15]*y[1] - x[16]*y[30] - x[17]*y[31] + x[18]*y[28] - x[19]*y[29] - x[20]*y[26] + x[21]*y[27] - x[22]*y[24] - x[23]*y[25] + x[24]*y[22] - x[25]*y[23] + x[26]*y[20] + x[27]*y[21] - x[28]*y[18] - x[29]*y[19] + x[30]*y[16] - x[31]*y[17],
        x[0]*y[15] + x[1]*y[14] - x[2]*y[13] + x[3]*y[12] + x[4]*y[11] - x[5]*y[10] + x[6]*y[9] + x[7]*y[8] - x[8]*y[7] + x[9]*y[6] - x[10]*y[5] - x[11]*y[4] + x[12]*y[3] + x[13]*y[2] - x[14]*y[1] + x[15]*y[0] + x[16]*y[31] - x[17]*y[30] + x[18]*y[29] + x[19]*y[28] - x[20]*y[27] - x[21]*y[26] + x[22]*y[25] - x[23]*y[24] + x[24]*y[23] + x[25]*y[22] - x[26]*y[21] + x[27]*y[20] + x[28]*y[19] - x[29]*y[18] + x[30]*y[17] + x[31]*y[16],
        x[0]*y[16] - x[1]*y[17] + x[2]*y[18] + x[3]*y[19] + x[4]*y[20] + x[5]*y[21] - x[6]*y[22] + x[7]*y[23] + x[8]*y[24] + x[9]*y[25] - x[10]*y[26] + x[11]*y[27] - x[12]*y[28] + x[13]*y[29] - x[14]*y[30] - x[15]*y[31] + x[16]*y[0] + x[17]*y[1] - x[18]*y[2] + x[19]*y[3] - x[20]*y[4] + x[21]*y[5] - x[22]*y[6] - x[23]*y[7] - x[24]*y[8] + x[25]*y[9] - x[26]*y[10] - x[27]*y[11] - x[28]*y[12] - x[29]*y[13] + x[30]*y[14] - x[31]*y[15],
        x[0]*y[17] + x[1]*y[16] - x[2]*y[19] + x[3]*y[18] - x[4]*y[21] + x[5]*y[20] - x[6]*y[23] - x[7]*y[22] - x[8]*y[25] + x[9]*y[24] - x[10]*y[27] - x[11]*y[26] - x[12]*y[29] - x[13]*y[28] + x[14]*y[31] - x[15]*y[30] - x[16]*y[1] + x[17]*y[0] - x[18]*y[3] - x[19]*y[2] - x[20]*y[5] - x[21]*y[4] + x[22]*y[7] - x[23]*y[6] - x[24]*y[9] - x[25]*y[8] + x[26]*y[11] - x[27]*y[10] + x[28]*y[13] - x[29]*y[12] + x[30]*y[15] + x[31]*y[14],
        x[0]*y[18] - x[1]*y[19] + x[2]*y[16] + x[3]*y[17] - x[4]*y[22] - x[5]*y[23] + x[6]*y[20] - x[7]*y[21] - x[8]*y[26] - x[9]*y[27] + x[10]*y[24] - x[11]*y[25] - x[12]*y[30] + x[13]*y[31] - x[14]*y[28] - x[15]*y[29] - x[16]*y[2] - x[17]*y[3] + x[18]*y[0] - x[19]*y[1] - x[20]*y[6] + x[21]*y[7] - x[22]*y[4] - x[23]*y[5] - x[24]*y[10] + x[25]*y[11] - x[26]*y[8] - x[27]*y[9] + x[28]*y[14] + x[29]*y[15] - x[30]*y[12] + x[31]*y[13],
        x[0]*y[19] + x[1]*y[18] - x[2]*y[17] + x[3]*y[16] + x[4]*y[23] - x[5]*y[22] + x[6]*y[21] + x[7]*y[20] + x[8]*y[27] - x[9]*y[26] + x[10]*y[25] + x[11]*y[24] - x[12]*y[31] - x[13]*y[30] + x[14]*y[29] - x[15]*y[28] + x[16]*y[3] - x[17]*y[2] + x[18]*y[1] + x[19]*y[0] - x[20]*y[7] - x[21]*y[6] + x[22]*y[5] - x[23]*y[4] - x[24]*y[11] - x[25]*y[10] + x[26]*y[9] - x[27]*y[8] - x[28]*y[15] + x[29]*y[14] - x[30]*y[13] - x[31]*y[12],
        x[0]*y[20] - x[1]*y[21] + x[2]*y[22] + x[3]*y[23] + x[4]*y[16] + x[5]*y[17] - x[6]*y[18] + x[7]*y[19] - x[8]*y[28] - x[9]*y[29] + x[10]*y[30] - x[11]*y[31] + x[12]*y[24] - x[13]*y[25] + x[14]*y[26] + x[15]*y[27] - x[16]*y[4] - x[17]*y[5] + x[18]*y[6] - x[19]*y[7] + x[20]*y[0] - x[21]*y[1] + x[22]*y[2] + x[23]*y[3] - x[24]*y[12] + x[25]*y[13] - x[26]*y[14] - x[27]*y[15] - x[28]*y[8] - x[29]*y[9] + x[30]*y[10] - x[31]*y[11],
        x[0]*y[21] + x[1]*y[20] - x[2]*y[23] + x[3]*y[22] - x[4]*y[17] + x[5]*y[16] - x[6]*y[19] - x[7]*y[18] + x[8]*y[29] - x[9]*y[28] + x[10]*y[31] + x[11]*y[30] + x[12]*y[25] + x[13]*y[24] - x[14]*y[27] + x[15]*y[26] + x[16]*y[5] - x[17]*y[4] + x[18]*y[7] + x[19]*y[6] + x[20]*y[1] + x[21]*y[0] - x[22]*y[3] + x[23]*y[2] - x[24]*y[13] - x[25]*y[12] + x[26]*y[15] - x[27]*y[14] + x[28]*y[9] - x[29]*y[8] + x[30]*y[11] + x[31]*y[10],
        x[0]*y[22] - x[1]*y[23] + x[2]*y[20] + x[3]*y[21] - x[4]*y[18] - x[5]*y[19] + x[6]*y[16] - x[7]*y[17] + x[8]*y[30] + x[9]*y[31] - x[10]*y[28] + x[11]*y[29] + x[12]*y[26] - x[13]*y[27] + x[14]*y[24] + x[15]*y[25] + x[16]*y[6] + x[17]*y[7] - x[18]*y[4] + x[19]*y[5] + x[20]*y[2] - x[21]*y[3] + x[22]*y[0] + x[23]*y[1] - x[24]*y[14] + x[25]*y[15] - x[26]*y[12] - x[27]*y[13] + x[28]*y[10] + x[29]*y[11] - x[30]*y[8] + x[31]*y[9],
        x[0]*y[23] + x[1]*y[22] - x[2]*y[21] + x[3]*y[20] + x[4]*y[19] - x[5]*y[18] + x[6]*y[17] + x[7]*y[16] - x[8]*y[31] + x[9]*y[30] - x[10]*y[29] - x[11]*y[28] + x[12]*y[27] + x[13]*y[26] - x[14]*y[25] + x[15]*y[24] - x[16]*y[7] + x[17]*y[6] - x[18]*y[5] - x[19]*y[4] + x[20]*y[3] + x[21]*y[2] - x[22]*y[1] + x[23]*y[0] - x[24]*y[15] - x[25]*y[14] + x[26]*y[13] - x[27]*y[12] - x[28]*y[11] + x[29]*y[10] - x[30]*y[9] - x[31]*y[8],
        x[0]*y[24] - x[1]*y[25] + x[2]*y[26] + x[3]*y[27] + x[4]*y[28] + x[5]*y[29] - x[6]*y[30] + x[7]*y[31] + x[8]*y[16] + x[9]*y[17] - x[10]*y[18] + x[11]*y[19] - x[12]*y[20] + x[13]*y[21] - x[14]*y[22] - x[15]*y[23] - x[16]*y[8] - x[17]*y[9] + x[18]*y[10] - x[19]*y[11] + x[20]*y[12] - x[21]*y[13] + x[22]*y[14] + x[23]*y[15] + x[24]*y[0] - x[25]*y[1] + x[26]*y[2] + x[27]*y[3] + x[28]*y[4] + x[29]*y[5] - x[30]*y[6] + x[31]*y[7],
        x[0]*y[25] + x[1]*y[24] - x[2]*y[27] + x[3]*y[26] - x[4]*y[29] + x[5]*y[28] - x[6]*y[31] - x[7]*y[30] - x[8]*y[17] + x[9]*y[16] - x[10]*y[19] - x[11]*y[18] - x[12]*y[21] - x[13]*y[20] + x[14]*y[23] - x[15]*y[22] + x[16]*y[9] - x[17]*y[8] + x[18]*y[11] + x[19]*y[10] + x[20]*y[13] + x[21]*y[12] - x[22]*y[15] + x[23]*y[14] + x[24]*y[1] + x[25]*y[0] - x[26]*y[3] + x[27]*y[2] - x[28]*y[5] + x[29]*y[4] - x[30]*y[7] - x[31]*y[6],
        x[0]*y[26] - x[1]*y[27] + x[2]*y[24] + x[3]*y[25] - x[4]*y[30] - x[5]*y[31] + x[6]*y[28] - x[7]*y[29] - x[8]*y[18] - x[9]*y[19] + x[10]*y[16] - x[11]*y[17] - x[12]*y[22] + x[13]*y[23] - x[14]*y[20] - x[15]*y[21] + x[16]*y[10] + x[17]*y[11] - x[18]*y[8] + x[19]*y[9] + x[20]*y[14] - x[21]*y[15] + x[22]*y[12] + x[23]*y[13] + x[24]*y[2] - x[25]*y[3] + x[26]*y[0] + x[27]*y[1] - x[28]*y[6] - x[29]*y[7] + x[30]*y[4] - x[31]*y[5],
        x[0]*y[27] + x[1]*y[26] - x[2]*y[25] + x[3]*y[24] + x[4]*y[31] - x[5]*y[30] + x[6]*y[29] + x[7]*y[28] + x[8]*y[19] - x[9]*y[18] + x[10]*y[17] + x[11]*y[16] - x[12]*y[23] - x[13]*y[22] + x[14]*y[21] - x[15]*y[20] - x[16]*y[11] + x[17]*y[10] - x[18]*y[9] - x[19]*y[8] + x[20]*y[15] + x[21]*y[14] - x[22]*y[13] + x[23]*y[12] + x[24]*y[3] + x[25]*y[2] - x[26]*y[1] + x[27]*y[0] + x[28]*y[7] - x[29]*y[6] + x[30]*y[5] + x[31]*y[4],
        x[0]*y[28] - x[1]*y[29] + x[2]*y[30] + x[3]*y[31] + x[4]*y[24] + x[5]*y[25] - x[6]*y[26] + x[7]*y[27] - x[8]*y[20] - x[9]*y[21] + x[10]*y[22] - x[11]*y[23] + x[12]*y[16] - x[13]*y[17] + x[14]*y[18] + x[15]*y[19] + x[16]*y[12] + x[17]*y[13] - x[18]*y[14] + x[19]*y[15] - x[20]*y[8] + x[21]*y[9] - x[22]*y[10] - x[23]*y[11] + x[24]*y[4] - x[25]*y[5] + x[26]*y[6] + x[27]*y[7] + x[28]*y[0] + x[29]*y[1] - x[30]*y[2] + x[31]*y[3],
        x[0]*y[29] + x[1]*y[28] - x[2]*y[31] + x[3]*y[30] - x[4]*y[25] + x[5]*y[24] - x[6]*y[27] - x[7]*y[26] + x[8]*y[21] - x[9]*y[20] + x[10]*y[23] + x[11]*y[22] + x[12]*y[17] + x[13]*y[16] - x[14]*y[19] + x[15]*y[18] - x[16]*y[13] + x[17]*y[12] - x[18]*y[15] - x[19]*y[14] - x[20]*y[9] - x[21]*y[8] + x[22]*y[11] - x[23]*y[10] + x[24]*y[5] + x[25]*y[4] - x[26]*y[7] + x[27]*y[6] - x[28]*y[1] + x[29]*y[0] - x[30]*y[3] - x[31]*y[2],
        x[0]*y[30] - x[1]*y[31] + x[2]*y[28] + x[3]*y[29] - x[4]*y[26] - x[5]*y[27] + x[6]*y[24] - x[7]*y[25] + x[8]*y[22] + x[9]*y[23] - x[10]*y[20] + x[11]*y[21] + x[12]*y[18] - x[13]*y[19] + x[14]*y[16] + x[15]*y[17] - x[16]*y[14] - x[17]*y[15] + x[18]*y[12] - x[19]*y[13] - x[20]*y[10] + x[21]*y[11] - x[22]*y[8] - x[23]*y[9] + x[24]*y[6] - x[25]*y[7] + x[26]*y[4] + x[27]*y[5] - x[28]*y[2] - x[29]*y[3] + x[30]*y[0] - x[31]*y[1],
        x[0]*y[31] + x[1]*y[30] - x[2]*y[29] + x[3]*y[28] + x[4]*y[27] - x[5]*y[26] + x[6]*y[25] + x[7]*y[24] - x[8]*y[23] + x[9]*y[22] - x[10]*y[21] - x[11]*y[20] + x[12]*y[19] + x[13]*y[18] - x[14]*y[17] + x[15]*y[16] + x[16]*y[15] - x[17]*y[14] + x[18]*y[13] + x[19]*y[12] - x[20]*y[11] - x[21]*y[10] + x[22]*y[9] - x[23]*y[8] + x[24]*y[7] + x[25]*y[6] - x[26]*y[5] + x[27]*y[4] + x[28]*y[3] - x[29]*y[2] + x[30]*y[1] + x[31]*y[0]});
}

static inline cga3
cga3_square(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){
        x[0]*x[0] - x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5] - x[6]*x[6] + x[7]*x[7] + x[8]*x[8] + x[9]*x[9] - x[10]*x[10] + x[11]*x[11] - x[12]*x[12] + x[13]*x[13] - x[14]*x[14] - x[15]*x[15] + x[16]*x[16] + x[17]*x[17] - x[18]*x[18] + x[19]*x[19] - x[20]*x[20] + x[21]*x[21] - x[22]*x[22] - x[23]*x[23] - x[24]*x[24] + x[25]*x[25] - x[26]*x[26] - x[27]*x[27] - x[28]*x[28] - x[29]*x[29] + x[30]*x[30] - x[31]*x[31],
        2*x[0]*x[1] - 2*x[6]*x[7] - 2*x[10]*x[11] - 2*x[12]*x[13] - 2*x[18]*x[19] - 2*x[20]*x[21] - 2*x[24]*x[25] + 2*x[30]*x[31],
        2*x[0]*x[2] - 2*x[5]*x[7] - 2*x[9]*x[11] - 2*x[12]*x[14] - 2*x[17]*x[19] - 2*x[20]*x[22] - 2*x[24]*x[26] + 2*x[29]*x[31],
        2*x[0]*x[3] + 2*x[4]*x[7] + 2*x[8]*x[11] - 2*x[12]*x[15] + 2*x[16]*x[19] - 2*x[20]*x[23] - 2*x[24]*x[27] - 2*x[28]*x[31],
        2*x[0]*x[4] + 2*x[3]*x[7] - 2*x[9]*x[13] + 2*x[10]*x[14] - 2*x[17]*x[21] + 2*x[18]*x[22] - 2*x[24]*x[28] - 2*x[27]*x[31],
        2*x[0]*x[5] - 2*x[2]*x[7] + 2*x[8]*x[13] + 2*x[10]*x[15] + 2*x[16]*x[21] + 2*x[18]*x[23] - 2*x[24]*x[29] + 2*x[26]*x[31],
        2*x[0]*x[6] - 2*x[1]*x[7] + 2*x[8]*x[14] + 2*x[9]*x[15] + 2*x[16]*x[22] + 2*x[17]*x[23] - 2*x[24]*x[30] + 2*x[25]*x[31],
        2*x[0]*x[7] + 2*x[1]*x[6] - 2*x[2]*x[5] + 2*x[3]*x[4] - 2*x[24]*x[31] - 2*x[25]*x[30] + 2*x[26]*x[29] - 2*x[27]*x[28],
        2*x[0]*x[8] + 2*x[3]*x[11] + 2*x[5]*x[13] - 2*x[6]*x[14] - 2*x[17]*x[25] + 2*x[18]*x[26] + 2*x[20]*x[28] + 2*x[23]*x[31],
        2*x[0]*x[9] - 2*x[2]*x[11] - 2*x[4]*x[13] - 2*x[6]*x[15] + 2*x[16]*x[25] + 2*x[18]*x[27] + 2*x[20]*x[29] - 2*x[22]*x[31],
        2*x[0]*x[10] - 2*x[1]*x[11] - 2*x[4]*x[14] - 2*x[5]*x[15] + 2*x[16]*x[26] + 2*x[17]*x[27] + 2*x[20]*x[30] - 2*x[21]*x[31],
        2*x[0]*x[11] + 2*x[1]*x[10] - 2*x[2]*x[9] + 2*x[3]*x[8] + 2*x[20]*x[31] + 2*x[21]*x[30] - 2*x[22]*x[29] + 2*x[23]*x[28],
        2*x[0]*x[12] - 2*x[1]*x[13] + 2*x[2]*x[14] + 2*x[3]*x[15] + 2*x[16]*x[28] + 2*x[17]*x[29] - 2*x[18]*x[30] + 2*x[19]*x[31],
        2*x[0]*x[13] + 2*x[1]*x[12] - 2*x[4]*x[9] + 2*x[5]*x[8] - 2*x[18]*x[31] - 2*x[19]*x[30] + 2*x[22]*x[27] - 2*x[23]*x[26],
        2*x[0]*x[14] + 2*x[2]*x[12] - 2*x[4]*x[10] + 2*x[6]*x[8] - 2*x[17]*x[31] - 2*x[19]*x[29] + 2*x[21]*x[27] - 2*x[23]*x[25],
        2*x[0]*x[15] + 2*x[3]*x[12] - 2*x[5]*x[10] + 2*x[6]*x[9] + 2*x[16]*x[31] + 2*x[19]*x[28] - 2*x[21]*x[26] + 2*x[22]*x[25],
        2*x[0]*x[16] + 2*x[3]*x[19] + 2*x[5]*x[21] - 2*x[6]*x[22] + 2*x[9]*x[25] - 2*x[10]*x[26] - 2*x[12]*x[28] - 2*x[15]*x[31],
        2*x[0]*x[17] - 2*x[2]*x[19] - 2*x[4]*x[21] - 2*x[6]*x[23] - 2*x[8]*x[25] - 2*x[10]*x[27] - 2*x[12]*x[29] + 2*x[14]*x[31],
        2*x[0]*x[18] - 2*x[1]*x[19] - 2*x[4]*x[22] - 2*x[5]*x[23] - 2*x[8]*x[26] - 2*x[9]*x[27] - 2*x[12]*x[30] + 2*x[13]*x[31],
        2*x[0]*x[19] + 2*x[1]*x[18] - 2*x[2]*x[17] + 2*x[3]*x[16] - 2*x[12]*x[31] - 2*x[13]*x[30] + 2*x[14]*x[29] - 2*x[15]*x[28],
        2*x[0]*x[20] - 2*x[1]*x[21] + 2*x[2]*x[22] + 2*x[3]*x[23] - 2*x[8]*x[28] - 2*x[9]*x[29] + 2*x[10]*x[30] - 2*x[11]*x[31],
        2*x[0]*x[21] + 2*x[1]*x[20] - 2*x[4]*x[17] + 2*x[5]*x[16] + 2*x[10]*x[31] + 2*x[11]*x[30] - 2*x[14]*x[27] + 2*x[15]*x[26],
        2*x[0]*x[22] + 2*x[2]*x[20] - 2*x[4]*x[18] + 2*x[6]*x[16] + 2*x[9]*x[31] + 2*x[11]*x[29] - 2*x[13]*x[27] + 2*x[15]*x[25],
        2*x[0]*x[23] + 2*x[3]*x[20] - 2*x[5]*x[18] + 2*x[6]*x[17] - 2*x[8]*x[31] - 2*x[11]*x[28] + 2*x[13]*x[26] - 2*x[14]*x[25],
        2*x[0]*x[24] - 2*x[1]*x[25] + 2*x[2]*x[26] + 2*x[3]*x[27] + 2*x[4]*x[28] + 2*x[5]*x[29] - 2*x[6]*x[30] + 2*x[7]*x[31],
        2*x[0]*x[25] + 2*x[1]*x[24] - 2*x[6]*x[31] - 2*x[7]*x[30] - 2*x[8]*x[17] + 2*x[9]*x[16] + 2*x[14]*x[23] - 2*x[15]*x[22],
        2*x[0]*x[26] + 2*x[2]*x[24] - 2*x[5]*x[31] - 2*x[7]*x[29] - 2*x[8]*x[18] + 2*x[10]*x[16] + 2*x[13]*x[23] - 2*x[15]*x[21],
        2*x[0]*x[27] + 2*x[3]*x[24] + 2*x[4]*x[31] + 2*x[7]*x[28] - 2*x[9]*x[18] + 2*x[10]*x[17] - 2*x[13]*x[22] + 2*x[14]*x[21],
        2*x[0]*x[28] + 2*x[3]*x[31] + 2*x[4]*x[24] + 2*x[7]*x[27] - 2*x[8]*x[20] - 2*x[11]*x[23] + 2*x[12]*x[16] + 2*x[15]*x[19],
        2*x[0]*x[29] - 2*x[2]*x[31] + 2*x[5]*x[24] - 2*x[7]*x[26] - 2*x[9]*x[20] + 2*x[11]*x[22] + 2*x[12]*x[17] - 2*x[14]*x[19],
        2*x[0]*x[30] - 2*x[1]*x[31] + 2*x[6]*x[24] - 2*x[7]*x[25] - 2*x[10]*x[20] + 2*x[11]*x[21] + 2*x[12]*x[18] - 2*x[13]*x[19],
        2*x[0]*x[31] + 2*x[1]*x[30] - 2*x[2]*x[29] + 2*x[3]*x[28] + 2*x[4]*x[27] - 2*x[5]*x[26] + 2*x[6]*x[25] + 2*x[7]*x[24] - 2*x[8]*x[23] + 2*x[9]*x[22] - 2*x[10]*x[21] - 2*x[11]*x[20] + 2*x[12]*x[19] + 2*x[13]*x[18] - 2*x[14]*x[17] + 2*x[15]*x[16]});
}

static inline cga3
cga3_conjugate(cga3 a)
{
    const double *x = a.data;
    return cga3_make((const double[]){x[0], -x[1], -x[2], -x[3], -x[4], -x[5], -x[6], x[7], -x[8], -x[9], -x[10], x[11], -x[12], x[13], x[14], x[15], -x[16], -x[17], -x[18], x[19], -x[20], x[21], x[22], x[23], -x[24], x[25], x[26], x[27], x[28], x[29], x[30], -x[31]});
}

static inline void
cga3_basis(cga3 result[32])
{
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) result[i].data[j] = i == j ? 1 : 0;
    }
}

static inline void
cga3_gens(cga3 result[5])
{
    result[0] = cga3_make((const double[ALGPY_CGA3_N]){[1] = 1});
    result[1] = cga3_make((const double[ALGPY_CGA3_N]){[2] = 1});
    result[2] = cga3_make((const double[ALGPY_CGA3_N]){[4] = 1});
    result[3] = cga3_make((const double[ALGPY_CGA3_N]){[8] = 1});
    result[4] = cga3_make((const double[ALGPY_CGA3_N]){[16] = 1});
}

static inline void
cga3_coeffs(cga3 result[1024])
{
    cga3 basis[32];
    cga3_basis(basis);
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) result[i*32 + j] = cga3_mul(basis[i], basis[j]);
    }
}

static inline double
cga3_scalar_product(cga3 a, cga3 b)
{
    const double *x = a.data, *y = b.data;
    return x[0]*y[0] - x[1]*y[1] + x[2]*y[2] + x[3]*y[3] + x[4]*y[4] + x[5]*y[5] - x[6]*y[6] + x[7]*y[7] + x[8]*y[8] + x[9]*y[9] - x[10]*y[10] + x[11]*y[11] - x[12]*y[12] + x[13]*y[13] - x[14]*y[14] - x[15]*y[15] + x[16]*y[16] + x[17]*y[17] - x[18]*y[18] + x[19]*y[19] - x[20]*y[20] + x[21]*y[21] - x[22]*y[22] - x[23]*y[23] - x[24]*y[24] + x[25]*y[25] - x[26]*y[26] - x[27]*y[27] - x[28]*y[28] - x[29]*y[29] + x[30]*y[30] - x[31]*y[31];
}

/* I = e0*e1*e2*e3*e4 is central and I*I = -1. */
static inline complex
cga3_center(cga3 a)
{
    return complex_make((const double[]){a.data[0], a.data[31]});
}

static inline cga3
cga3_from_center(complex z)
{
    cga3 result = {{0}};
    result.data[0] = z.data[0];
    result.data[31] = z.data[1];
    return result;
}

static inline cga3
cga3_center_scale(cga3 a, complex z)
{
    const double *x = a.data;
    double re = z.data[0], im = z.data[1];
    return cga3_make((const double[]){re*x[0] - im*x[31], re*x[1] + im*x[30], re*x[2] + im*x[29], re*x[3] - im*x[28], re*x[4] - im*x[27], re*x[5] + im*x[26], re*x[6] + im*x[25], re*x[7] - im*x[24], re*x[8] + im*x[23], re*x[9] - im*x[22], re*x[10] - im*x[21], re*x[11] + im*x[20], re*x[12] + im*x[19], re*x[13] - im*x[18], re*x[14] - im*x[17], re*x[15] + im*x[16], re*x[16] - im*x[15], re*x[17] + im*x[14], re*x[18] + im*x[13], re*x[19] - im*x[12], re*x[20] - im*x[11], re*x[21] + im*x[10], re*x[22] + im*x[9], re*x[23] - im*x[8], re*x[24] + im*x[7], re*x[25] - im*x[6], re*x[26] - im*x[5], re*x[27] + im*x[4], re*x[28] + im*x[3], re*x[29] - im*x[2], re*x[30] - im*x[1], re*x[31] + im*x[0]});
}

static inline complex
cga3_center_product(cga3 a, cga3 b)
{
    const double *x = a.data, *y = b.data;
    return complex_make((const double[]){cga3_scalar_product(a, b),
        x[0]*y[31] + x[1]*y[30] - x[2]*y[29] + x[3]*y[28] + x[4]*y[27] - x[5]*y[26] + x[6]*y[25] + x[7]*y[24] - x[8]*y[23] + x[9]*y[22] - x[10]*y[21] - x[11]*y[20] + x[12]*y[19] + x[13]*y[18] - x[14]*y[17] + x[15]*y[16] + x[16]*y[15] - x[17]*y[14] + x[18]*y[13] + x[19]*y[12] - x[20]*y[11] - x[21]*y[10] + x[22]*y[9] - x[23]*y[8] + x[24]*y[7] + x[25]*y[6] - x[26]*y[5] + x[27]*y[4] + x[28]*y[3] - x[29]*y[2] + x[30]*y[1] + x[31]*y[0]});
}

/* The first 16 blades form a basis over the complex centre. */
static inline complex
cga3_coordinate(cga3 a, int i)
{
    static const double signs[16] = {1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1};
    return complex_make((const double[]){a.data[i], signs[i]*a.data[31 ^ i]});
}

/* Normalized complex traces give a quartic, rather than a real octic. */
static inline void
cga3_centered_quartic(cga3 b2, cga3 b3, complex p[5])
{
    complex s2 = cga3_center(b2), s3 = cga3_center(b3);
    complex s4 = cga3_center_product(b2, b2);
    p[0] = vec2_sub(vec2_scale(complex_square(s2), 2), s4);
    p[1] = vec2_scale(s3, -4./3);
    p[2] = vec2_scale(s2, -2);
    p[3] = complex_from_re(0);
    p[4] = complex_from_re(1);
}

static inline void
cga3_complex_minpoly(cga3 a, complex p[5])
{
    complex w = cga3_center(a), w2 = complex_square(w), q[5];
    cga3 b = a;
    b.data[0] = b.data[31] = 0;
    cga3 b2 = cga3_square(b), b3 = cga3_mul(b2, b);
    cga3_centered_quartic(b2, b3, q);
    p[0] = vec2_add(vec2_sub(complex_square(w2), complex_mul(w, q[1])),
        vec2_add(complex_mul(w2, q[2]), q[0]));
    p[1] = vec2_sub(vec2_sub(q[1], complex_mul(vec2_scale(w, 2), q[2])), complex_mul(vec2_scale(w, 4), w2));
    p[2] = vec2_add(vec2_scale(w2, 6), q[2]);
    p[3] = vec2_scale(w, -4);
    p[4] = complex_from_re(1);
}

/* Generic REAL minimal polynomial: p(x)*conj(p)(x), degree eight. */
static inline void
cga3_minpoly(cga3 a, double result[9])
{
    complex p[5];
    cga3_complex_minpoly(a, p);
    for (int i = 0; i < 9; i++) result[i] = 0;
    for (int i = 0; i < 5; i++) {
        result[2*i] += p[i].data[0]*p[i].data[0] + p[i].data[1]*p[i].data[1];
        for (int j = 0; j < i; j++) {
            result[i + j] += 2*(p[i].data[0]*p[j].data[0] + p[i].data[1]*p[j].data[1]);
        }
    }
}

static inline double
cga3_norm8(cga3 a)
{
    complex p[5];
    cga3_complex_minpoly(a, p);
    return complex_norm2(p[0]);
}

static inline double
cga3_norm(cga3 a)
{
    complex p[5];
    cga3_complex_minpoly(a, p);
    return sqrt(sqrt(hypot(p[0].data[0], p[0].data[1])));
}

static inline cga3
cga3_inv(cga3 a)
{
    complex p[5];
    cga3_complex_minpoly(a, p);
    cga3 a2 = cga3_square(a), a3 = cga3_mul(a2, a);
    cga3 numerator = vec32_add(a3, vec32_add(cga3_center_scale(a2, p[3]),
        vec32_add(cga3_center_scale(a, p[2]), cga3_from_center(p[1]))));
    return cga3_center_scale(numerator, vec2_neg(complex_inv(p[0])));
}

static inline void
cga3_eivals(cga3 a, complex result[8])
{
    complex w = cga3_center(a), p[5];
    cga3 b = a;
    b.data[0] = b.data[31] = 0;
    cga3 b2 = cga3_square(b), b3 = cga3_mul(b2, b);
    cga3_centered_quartic(b2, b3, p);
    poly4_complex_coeffs_to_roots(p, result);
    for (int i = 0; i < 4; i++) {
        result[i] = vec2_add(w, result[i]);
        result[4 + i] = complex_conjugate(result[i]);
    }
}

static inline cga3
cga3_div(cga3 a, cga3 b)
{
    return cga3_mul(a, cga3_inv(b));
}

static inline cga3
cga3_nan(void)
{
    cga3 result;
    for (int i = 0; i < 32; i++) result.data[i] = NAN;
    return result;
}

static inline int
cga3_unused_generator(cga3 a)
{
    for (int generator = 1; generator < 5; generator++) {
        int used = 0;
        for (int i = 1; i < 32; i++) {
            if ((i & (1 << generator)) && a.data[i] != 0) { used = 1; break; }
        }
        if (!used) return generator;
    }
    return -1;
}

static inline cga3
cga3_apply_subalgebra(cga3 a, int generator,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex),
    complex (*cderiv2)(complex), complex (*cderiv3)(complex))
{
    cga2 suba;
    cga3 result = {{0}};
    int low = (1 << generator) - 1;
    (void)cderiv2; (void)cderiv3;
    for (int i = 0; i < 16; i++) suba.data[i] = a.data[(i & low) | ((i & ~low) << 1)];
    cga2 value = cga2_apply(suba, f, deriv1, deriv2, deriv3, cf, cderiv1);
    if (isnan(value.data[0])) return cga3_nan();
    for (int i = 0; i < 16; i++) result.data[(i & low) | ((i & ~low) << 1)] = value.data[i];
    return result;
}

static inline int
cga3_quadratic_relation(cga3 b, cga3 b2, complex *linear)
{
    int pivot = 1;
    double size = 0;
    for (int i = 2; i < 16; i++) {
        if (poly_complex_abs(cga3_coordinate(b, i)) > poly_complex_abs(cga3_coordinate(b, pivot))) pivot = i;
    }
    *linear = complex_div(cga3_coordinate(b2, pivot), cga3_coordinate(b, pivot));
    if (!hermite_complex_is_finite(*linear)) return 0;
    cga3 product = cga3_center_scale(b, *linear);
    for (int i = 0; i < 32; i++) size = fmax(size, fabs(b2.data[i]) + fabs(product.data[i]));
    for (int i = 1; i < 31; i++) {
        if (fabs(b2.data[i] - product.data[i]) > 32*DBL_EPSILON*size) return 0;
    }
    return 1;
}

static inline int
cga3_cubic_relation(cga3 b, cga3 b2, cga3 b3, complex p[4])
{
    int first = 1, second = 1;
    double largest = 0, bound = 0;
    for (int i = 2; i < 16; i++) {
        if (poly_complex_abs(cga3_coordinate(b, i)) > poly_complex_abs(cga3_coordinate(b, first))) first = i;
    }
    complex ratio = complex_div(cga3_coordinate(b2, first), cga3_coordinate(b, first));
    for (int i = 1; i < 16; i++) {
        double residual = poly_complex_abs(vec2_sub(cga3_coordinate(b2, i), complex_mul(ratio, cga3_coordinate(b, i))));
        if (residual > largest) { largest = residual; second = i; }
    }
    if (largest == 0) return 0;
    complex quadratic = complex_div(vec2_sub(cga3_coordinate(b3, second),
        complex_mul(complex_div(cga3_coordinate(b3, first), cga3_coordinate(b, first)), cga3_coordinate(b, second))),
        vec2_sub(cga3_coordinate(b2, second), complex_mul(ratio, cga3_coordinate(b, second))));
    complex linear = complex_div(vec2_sub(cga3_coordinate(b3, first),
        complex_mul(quadratic, cga3_coordinate(b2, first))), cga3_coordinate(b, first));
    if (!hermite_complex_is_finite(linear) || !hermite_complex_is_finite(quadratic)) return 0;
    cga3 q = cga3_center_scale(b2, quadratic), l = cga3_center_scale(b, linear);
    for (int i = 0; i < 32; i++) bound = fmax(bound, fabs(b3.data[i]) + fabs(q.data[i]) + fabs(l.data[i]));
    for (int i = 1; i < 31; i++) {
        if (fabs(b3.data[i] - q.data[i] - l.data[i]) > 64*DBL_EPSILON*bound) return 0;
    }
    p[0] = vec2_sub(complex_mul(quadratic, cga3_center(b2)), cga3_center(b3));
    p[1] = vec2_neg(linear);
    p[2] = vec2_neg(quadratic);
    p[3] = complex_from_re(1);
    return 1;
}

/* Quartic calculus over span{1,I}; no general Clifford dependency. */
static inline cga3
cga3_apply(cga3 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex),
    complex (*cderiv2)(complex), complex (*cderiv3)(complex))
{
    complex w = cga3_center(a), linear, coefficients[4], p[5];
    cga3 b = a;
    b.data[0] = b.data[31] = 0;
    if (vec32_is_real(a)) {
        complex value = cf(w);
        cga3 result = cga3_from_re(value.data[0]);
        result.data[1] = value.data[1];
        return result;
    }
    if (!vec32_any_nonzero(b)) return cga3_from_center(cf(w));
    int generator = cga3_unused_generator(a);
    if (generator >= 0) {
        return cga3_apply_subalgebra(a, generator,
            f, deriv1, deriv2, deriv3, cf, cderiv1, cderiv2, cderiv3);
    }
    cga3 b2 = cga3_square(b), b3;
    int ok;
    if (cga3_quadratic_relation(b, b2, &linear)) {
        complex shift = vec2_scale(linear, .5);
        complex root = complex_sqrt(vec2_add(cga3_center(b2), complex_square(shift)));
        complex nodes[2] = {vec2_sub(shift, root), vec2_add(shift, root)};
        ok = hermite_complex_interpolate(2, w, nodes,
            f, cf, cderiv1, cderiv2, cderiv3, coefficients);
        if (!ok) return cga3_nan();
        return vec32_add(cga3_from_center(coefficients[0]), cga3_center_scale(b, coefficients[1]));
    }
    b3 = cga3_mul(b2, b);
    if (cga3_cubic_relation(b, b2, b3, p)) {
        complex roots[3];
        poly3_complex_coeffs_to_roots(p, roots);
        ok = hermite_complex_interpolate(3, w, roots,
            f, cf, cderiv1, cderiv2, cderiv3, coefficients);
    } else {
        cga3_centered_quartic(b2, b3, p);
        ok = poly4_complex_apply(p, w, f, cf, cderiv1, cderiv2, cderiv3, coefficients);
    }
    if (!ok) return cga3_nan();
    return vec32_add(cga3_from_center(coefficients[0]),
        vec32_add(cga3_center_scale(b, coefficients[1]),
            vec32_add(cga3_center_scale(b2, coefficients[2]), cga3_center_scale(b3, coefficients[3]))));
}

static inline cga3
cga3_exp(cga3 a)
{
    return cga3_apply(a, double_exp, double_exp_deriv1,
        double_exp_deriv2, double_exp_deriv3, complex_exp, complex_exp_deriv1,
        complex_exp_deriv2, complex_exp_deriv3);
}

static inline cga3
cga3_expm1(cga3 a)
{
    return cga3_apply(a, double_expm1, double_expm1_deriv1,
        double_expm1_deriv2, double_expm1_deriv3, complex_expm1, complex_expm1_deriv1,
        complex_expm1_deriv2, complex_expm1_deriv3);
}

static inline cga3
cga3_log(cga3 a)
{
    return cga3_apply(a, double_log, double_log_deriv1,
        double_log_deriv2, double_log_deriv3, complex_log, complex_log_deriv1,
        complex_log_deriv2, complex_log_deriv3);
}

static inline cga3
cga3_log1p(cga3 a)
{
    return cga3_apply(a, double_log1p, double_log1p_deriv1,
        double_log1p_deriv2, double_log1p_deriv3, complex_log1p, complex_log1p_deriv1,
        complex_log1p_deriv2, complex_log1p_deriv3);
}

static inline cga3
cga3_sqrt(cga3 a)
{
    return cga3_apply(a, double_sqrt, double_sqrt_deriv1,
        double_sqrt_deriv2, double_sqrt_deriv3, complex_sqrt, complex_sqrt_deriv1,
        complex_sqrt_deriv2, complex_sqrt_deriv3);
}

static inline cga3
cga3_sin(cga3 a)
{
    return cga3_apply(a, double_sin, double_sin_deriv1,
        double_sin_deriv2, double_sin_deriv3, complex_sin, complex_sin_deriv1,
        complex_sin_deriv2, complex_sin_deriv3);
}

static inline cga3
cga3_cos(cga3 a)
{
    return cga3_apply(a, double_cos, double_cos_deriv1,
        double_cos_deriv2, double_cos_deriv3, complex_cos, complex_cos_deriv1,
        complex_cos_deriv2, complex_cos_deriv3);
}

static inline cga3
cga3_tan(cga3 a)
{
    return cga3_apply(a, double_tan, double_tan_deriv1,
        double_tan_deriv2, double_tan_deriv3, complex_tan, complex_tan_deriv1,
        complex_tan_deriv2, complex_tan_deriv3);
}

static inline cga3
cga3_asin(cga3 a)
{
    return cga3_apply(a, double_asin, double_asin_deriv1,
        double_asin_deriv2, double_asin_deriv3, complex_asin, complex_asin_deriv1,
        complex_asin_deriv2, complex_asin_deriv3);
}

static inline cga3
cga3_acos(cga3 a)
{
    return cga3_apply(a, double_acos, double_acos_deriv1,
        double_acos_deriv2, double_acos_deriv3, complex_acos, complex_acos_deriv1,
        complex_acos_deriv2, complex_acos_deriv3);
}

static inline cga3
cga3_atan(cga3 a)
{
    return cga3_apply(a, double_atan, double_atan_deriv1,
        double_atan_deriv2, double_atan_deriv3, complex_atan, complex_atan_deriv1,
        complex_atan_deriv2, complex_atan_deriv3);
}

static inline cga3
cga3_sinh(cga3 a)
{
    return cga3_apply(a, double_sinh, double_sinh_deriv1,
        double_sinh_deriv2, double_sinh_deriv3, complex_sinh, complex_sinh_deriv1,
        complex_sinh_deriv2, complex_sinh_deriv3);
}

static inline cga3
cga3_cosh(cga3 a)
{
    return cga3_apply(a, double_cosh, double_cosh_deriv1,
        double_cosh_deriv2, double_cosh_deriv3, complex_cosh, complex_cosh_deriv1,
        complex_cosh_deriv2, complex_cosh_deriv3);
}

static inline cga3
cga3_tanh(cga3 a)
{
    return cga3_apply(a, double_tanh, double_tanh_deriv1,
        double_tanh_deriv2, double_tanh_deriv3, complex_tanh, complex_tanh_deriv1,
        complex_tanh_deriv2, complex_tanh_deriv3);
}

static inline cga3
cga3_asinh(cga3 a)
{
    return cga3_apply(a, double_asinh, double_asinh_deriv1,
        double_asinh_deriv2, double_asinh_deriv3, complex_asinh, complex_asinh_deriv1,
        complex_asinh_deriv2, complex_asinh_deriv3);
}

static inline cga3
cga3_acosh(cga3 a)
{
    return cga3_apply(a, double_acosh, double_acosh_deriv1,
        double_acosh_deriv2, double_acosh_deriv3, complex_acosh, complex_acosh_deriv1,
        complex_acosh_deriv2, complex_acosh_deriv3);
}

static inline cga3
cga3_atanh(cga3 a)
{
    return cga3_apply(a, double_atanh, double_atanh_deriv1,
        double_atanh_deriv2, double_atanh_deriv3, complex_atanh, complex_atanh_deriv1,
        complex_atanh_deriv2, complex_atanh_deriv3);
}

static inline cga3
cga3_radd(double a, cga3 b)
{
    return cga3_add(cga3_from_re(a), b);
}

static inline cga3
cga3_ladd(cga3 a, double b)
{
    return cga3_add(a, cga3_from_re(b));
}

static inline cga3
cga3_rsub(double a, cga3 b)
{
    return cga3_sub(cga3_from_re(a), b);
}

static inline cga3
cga3_lsub(cga3 a, double b)
{
    return cga3_sub(a, cga3_from_re(b));
}

static inline cga3
cga3_rmul(double a, cga3 b)
{
    return cga3_mul(cga3_from_re(a), b);
}

static inline cga3
cga3_lmul(cga3 a, double b)
{
    return cga3_mul(a, cga3_from_re(b));
}

static inline cga3
cga3_rdiv(double a, cga3 b)
{
    return cga3_div(cga3_from_re(a), b);
}

static inline cga3
cga3_ldiv(cga3 a, double b)
{
    return cga3_div(a, cga3_from_re(b));
}

static inline cga3
cga3_ipow(cga3 a, int64_t b)
{
    cga3 result = cga3_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return cga3_mul(a, a);
    }
    if (b == -1) {
        return cga3_inv(a);
    }
    if (b < 0) {
        a = cga3_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = cga3_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = cga3_mul(a, a);
        }
    }
    return result;
}

static inline cga3
cga3_rpow(cga3 a, cga3 b)
{
    return cga3_exp(cga3_mul(cga3_log(a), b));
}

static inline cga3
cga3_pow(cga3 a, cga3 b)
{
    return cga3_rpow(a, b);
}

static inline cga3
cga3_lpow(cga3 a, cga3 b)
{
    return cga3_exp(cga3_mul(b, cga3_log(a)));
}

static inline cga3
cga3_rpows(double a, cga3 b)
{
    return cga3_rpow(cga3_from_re(a), b);
}

static inline cga3
cga3_lpows(cga3 a, double b)
{
    return cga3_rpow(a, cga3_from_re(b));
}

static inline cga3
cga3_mix(cga3 a, cga3 b, double t)
{
    cga3 delta = cga3_mul(cga3_inv(a), b);
    return cga3_mul(a, cga3_lpows(delta, t));
}

static inline void
cga3_lmat(cga3 a, double result[ALGPY_CGA3_N*ALGPY_CGA3_N])
{
    cga3 basis = {{0}};
    cga3 product;
    int row, col;
    for (col = 0; col < ALGPY_CGA3_N; col++) {
        basis.data[col] = 1;
        product = cga3_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CGA3_N; row++) {
            result[row * ALGPY_CGA3_N + col] = product.data[row];
        }
    }
}

static inline void
cga3_rmat(cga3 a, double result[ALGPY_CGA3_N*ALGPY_CGA3_N])
{
    cga3 basis = {{0}};
    cga3 product;
    int row, col;
    for (row = 0; row < ALGPY_CGA3_N; row++) {
        basis.data[row] = 1;
        product = cga3_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CGA3_N; col++) {
            result[row * ALGPY_CGA3_N + col] = product.data[col];
        }
    }
}

/* Lattice API. The exterior and regressive products and the Poincare duals
 * do not depend on the metric, so they are those of exterior5, which shares
 * this coefficient order: join with identity bottom = 1, meet with identity
 * top = the pseudoscalar, rdual(e_a) join e_a = top = e_a join ldual(e_a). */
static inline cga3
cga3_bottom(void)
{
    return exterior5_bottom();
}

static inline cga3
cga3_top(void)
{
    return exterior5_top();
}

static inline cga3
cga3_join(cga3 a, cga3 b)
{
    return exterior5_join(a, b);
}

static inline cga3
cga3_meet(cga3 a, cga3 b)
{
    return exterior5_meet(a, b);
}

static inline cga3
cga3_ldual(cga3 a)
{
    return exterior5_ldual(a);
}

static inline cga3
cga3_rdual(cga3 a)
{
    return exterior5_rdual(a);
}

static inline cga3
cga3_dual(cga3 a)
{
    return exterior5_dual(a);
}

/*
 * Minimal polynomial over the centre span{1, w}, w = e01234 with w*w = -1,
 * which is the complex1 algebra. Coefficients are in ascending order, p[i]
 * standing for p[i][0] + p[i][1]*w.
 */
static inline void
cga3_center_minpoly(cga3 a, complex1 p[5])
{
    cga3_complex_minpoly(a, p);
}

#endif /* ALGPY_CGA3_H */
