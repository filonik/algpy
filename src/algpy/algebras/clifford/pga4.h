/* Projective geometric algebra in 4D.
 * Generator squares: 0, 1, 1, 1, 1.
 * 5 generators; 32 real coefficients. */
#ifndef ALGPY_PGA4_H
#define ALGPY_PGA4_H

#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../../polynomials/hermite8.h"
#include "../../vectorspaces/vec32.h"
#include "../../lattices/clifford/exterior5.h"
#include "../multicomplex/dual1.h"
#include "pga3.h"

#define ALGPY_PGA4_N 32

typedef vec32 pga4;

static inline pga4
pga4_make(const double coeffs[ALGPY_PGA4_N])
{
    pga4 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline pga4
pga4_grade0(pga4 a)
{
    return pga4_make((const double[ALGPY_PGA4_N]){[0] = a.data[0]});
}

static inline pga4
pga4_grade1(pga4 a)
{
    return pga4_make((const double[ALGPY_PGA4_N]){[1] = a.data[1], [2] = a.data[2], [4] = a.data[4], [8] = a.data[8], [16] = a.data[16]});
}

static inline pga4
pga4_grade2(pga4 a)
{
    return pga4_make((const double[ALGPY_PGA4_N]){[3] = a.data[3], [5] = a.data[5], [6] = a.data[6], [9] = a.data[9], [10] = a.data[10], [12] = a.data[12], [17] = a.data[17], [18] = a.data[18], [20] = a.data[20], [24] = a.data[24]});
}

static inline pga4
pga4_grade3(pga4 a)
{
    return pga4_make((const double[ALGPY_PGA4_N]){[7] = a.data[7], [11] = a.data[11], [13] = a.data[13], [14] = a.data[14], [19] = a.data[19], [21] = a.data[21], [22] = a.data[22], [25] = a.data[25], [26] = a.data[26], [28] = a.data[28]});
}

static inline pga4
pga4_grade4(pga4 a)
{
    return pga4_make((const double[ALGPY_PGA4_N]){[15] = a.data[15], [23] = a.data[23], [27] = a.data[27], [29] = a.data[29], [30] = a.data[30]});
}

static inline pga4
pga4_grade5(pga4 a)
{
    return pga4_make((const double[ALGPY_PGA4_N]){[31] = a.data[31]});
}

static inline pga4
pga4_from_re(double a)
{
    return pga4_make((const double[ALGPY_PGA4_N]){a});
}

static inline double
pga4_to_re(pga4 a)
{
    return a.data[0];
}

static inline pga4
pga4_zero(void)
{
    return vec32_zero();
}

static inline pga4
pga4_one(void)
{
    return pga4_from_re(1);
}

static inline pga4
pga4_pos(pga4 a)
{
    return vec32_pos(a);
}

static inline pga4
pga4_neg(pga4 a)
{
    return vec32_neg(a);
}

static inline pga4
pga4_add(pga4 a, pga4 b)
{
    return vec32_add(a, b);
}

static inline pga4
pga4_sub(pga4 a, pga4 b)
{
    return vec32_sub(a, b);
}

static inline int
pga4_equal(pga4 a, pga4 b)
{
    return vec32_all_eq(a, b);
}

static inline int
pga4_not_equal(pga4 a, pga4 b)
{
    return !pga4_equal(a, b);
}

static inline int
pga4_is_real(pga4 a)
{
    return vec32_is_real(a);
}

static inline pga4
pga4_re(pga4 a)
{
    return pga4_from_re(a.data[0]);
}

static inline pga4
pga4_im(pga4 a)
{
    return pga4_sub(a, pga4_re(a));
}

static inline pga4
pga4_mul(pga4 a, pga4 b)
{
    if (vec32_is_real(b)) return vec32_scale(a, b.data[0]);
    if (vec32_is_real(a)) return vec32_scale(b, a.data[0]);
    const double *x = a.data, *y = b.data;
    return pga4_make((const double[]){
        x[0]*y[0] + x[2]*y[2] + x[4]*y[4] - x[6]*y[6] + x[8]*y[8] - x[10]*y[10] - x[12]*y[12] - x[14]*y[14] + x[16]*y[16] - x[18]*y[18] - x[20]*y[20] - x[22]*y[22] - x[24]*y[24] - x[26]*y[26] - x[28]*y[28] + x[30]*y[30],
        x[0]*y[1] + x[1]*y[0] - x[2]*y[3] + x[3]*y[2] - x[4]*y[5] + x[5]*y[4] - x[6]*y[7] - x[7]*y[6] - x[8]*y[9] + x[9]*y[8] - x[10]*y[11] - x[11]*y[10] - x[12]*y[13] - x[13]*y[12] + x[14]*y[15] - x[15]*y[14] - x[16]*y[17] + x[17]*y[16] - x[18]*y[19] - x[19]*y[18] - x[20]*y[21] - x[21]*y[20] + x[22]*y[23] - x[23]*y[22] - x[24]*y[25] - x[25]*y[24] + x[26]*y[27] - x[27]*y[26] + x[28]*y[29] - x[29]*y[28] + x[30]*y[31] + x[31]*y[30],
        x[0]*y[2] + x[2]*y[0] - x[4]*y[6] + x[6]*y[4] - x[8]*y[10] + x[10]*y[8] - x[12]*y[14] - x[14]*y[12] - x[16]*y[18] + x[18]*y[16] - x[20]*y[22] - x[22]*y[20] - x[24]*y[26] - x[26]*y[24] + x[28]*y[30] - x[30]*y[28],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0] + x[4]*y[7] - x[5]*y[6] + x[6]*y[5] + x[7]*y[4] + x[8]*y[11] - x[9]*y[10] + x[10]*y[9] + x[11]*y[8] - x[12]*y[15] - x[13]*y[14] + x[14]*y[13] - x[15]*y[12] + x[16]*y[19] - x[17]*y[18] + x[18]*y[17] + x[19]*y[16] - x[20]*y[23] - x[21]*y[22] + x[22]*y[21] - x[23]*y[20] - x[24]*y[27] - x[25]*y[26] + x[26]*y[25] - x[27]*y[24] - x[28]*y[31] + x[29]*y[30] - x[30]*y[29] - x[31]*y[28],
        x[0]*y[4] + x[2]*y[6] + x[4]*y[0] - x[6]*y[2] - x[8]*y[12] + x[10]*y[14] + x[12]*y[8] + x[14]*y[10] - x[16]*y[20] + x[18]*y[22] + x[20]*y[16] + x[22]*y[18] - x[24]*y[28] - x[26]*y[30] - x[28]*y[24] + x[30]*y[26],
        x[0]*y[5] + x[1]*y[4] - x[2]*y[7] + x[3]*y[6] - x[4]*y[1] + x[5]*y[0] - x[6]*y[3] - x[7]*y[2] + x[8]*y[13] - x[9]*y[12] + x[10]*y[15] + x[11]*y[14] + x[12]*y[9] + x[13]*y[8] - x[14]*y[11] + x[15]*y[10] + x[16]*y[21] - x[17]*y[20] + x[18]*y[23] + x[19]*y[22] + x[20]*y[17] + x[21]*y[16] - x[22]*y[19] + x[23]*y[18] - x[24]*y[29] - x[25]*y[28] + x[26]*y[31] - x[27]*y[30] + x[28]*y[25] - x[29]*y[24] + x[30]*y[27] + x[31]*y[26],
        x[0]*y[6] + x[2]*y[4] - x[4]*y[2] + x[6]*y[0] + x[8]*y[14] - x[10]*y[12] + x[12]*y[10] + x[14]*y[8] + x[16]*y[22] - x[18]*y[20] + x[20]*y[18] + x[22]*y[16] - x[24]*y[30] - x[26]*y[28] + x[28]*y[26] - x[30]*y[24],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0] - x[8]*y[15] + x[9]*y[14] - x[10]*y[13] - x[11]*y[12] + x[12]*y[11] + x[13]*y[10] - x[14]*y[9] + x[15]*y[8] - x[16]*y[23] + x[17]*y[22] - x[18]*y[21] - x[19]*y[20] + x[20]*y[19] + x[21]*y[18] - x[22]*y[17] + x[23]*y[16] - x[24]*y[31] - x[25]*y[30] + x[26]*y[29] - x[27]*y[28] - x[28]*y[27] + x[29]*y[26] - x[30]*y[25] - x[31]*y[24],
        x[0]*y[8] + x[2]*y[10] + x[4]*y[12] - x[6]*y[14] + x[8]*y[0] - x[10]*y[2] - x[12]*y[4] - x[14]*y[6] - x[16]*y[24] + x[18]*y[26] + x[20]*y[28] + x[22]*y[30] + x[24]*y[16] + x[26]*y[18] + x[28]*y[20] - x[30]*y[22],
        x[0]*y[9] + x[1]*y[8] - x[2]*y[11] + x[3]*y[10] - x[4]*y[13] + x[5]*y[12] - x[6]*y[15] - x[7]*y[14] - x[8]*y[1] + x[9]*y[0] - x[10]*y[3] - x[11]*y[2] - x[12]*y[5] - x[13]*y[4] + x[14]*y[7] - x[15]*y[6] + x[16]*y[25] - x[17]*y[24] + x[18]*y[27] + x[19]*y[26] + x[20]*y[29] + x[21]*y[28] - x[22]*y[31] + x[23]*y[30] + x[24]*y[17] + x[25]*y[16] - x[26]*y[19] + x[27]*y[18] - x[28]*y[21] + x[29]*y[20] - x[30]*y[23] - x[31]*y[22],
        x[0]*y[10] + x[2]*y[8] - x[4]*y[14] + x[6]*y[12] - x[8]*y[2] + x[10]*y[0] - x[12]*y[6] - x[14]*y[4] + x[16]*y[26] - x[18]*y[24] + x[20]*y[30] + x[22]*y[28] + x[24]*y[18] + x[26]*y[16] - x[28]*y[22] + x[30]*y[20],
        x[0]*y[11] + x[1]*y[10] - x[2]*y[9] + x[3]*y[8] + x[4]*y[15] - x[5]*y[14] + x[6]*y[13] + x[7]*y[12] + x[8]*y[3] - x[9]*y[2] + x[10]*y[1] + x[11]*y[0] - x[12]*y[7] - x[13]*y[6] + x[14]*y[5] - x[15]*y[4] - x[16]*y[27] + x[17]*y[26] - x[18]*y[25] - x[19]*y[24] + x[20]*y[31] + x[21]*y[30] - x[22]*y[29] + x[23]*y[28] + x[24]*y[19] + x[25]*y[18] - x[26]*y[17] + x[27]*y[16] + x[28]*y[23] - x[29]*y[22] + x[30]*y[21] + x[31]*y[20],
        x[0]*y[12] + x[2]*y[14] + x[4]*y[8] - x[6]*y[10] - x[8]*y[4] + x[10]*y[6] + x[12]*y[0] + x[14]*y[2] + x[16]*y[28] - x[18]*y[30] - x[20]*y[24] - x[22]*y[26] + x[24]*y[20] + x[26]*y[22] + x[28]*y[16] - x[30]*y[18],
        x[0]*y[13] + x[1]*y[12] - x[2]*y[15] + x[3]*y[14] - x[4]*y[9] + x[5]*y[8] - x[6]*y[11] - x[7]*y[10] + x[8]*y[5] - x[9]*y[4] + x[10]*y[7] + x[11]*y[6] + x[12]*y[1] + x[13]*y[0] - x[14]*y[3] + x[15]*y[2] - x[16]*y[29] + x[17]*y[28] - x[18]*y[31] - x[19]*y[30] - x[20]*y[25] - x[21]*y[24] + x[22]*y[27] - x[23]*y[26] + x[24]*y[21] + x[25]*y[20] - x[26]*y[23] + x[27]*y[22] - x[28]*y[17] + x[29]*y[16] - x[30]*y[19] - x[31]*y[18],
        x[0]*y[14] + x[2]*y[12] - x[4]*y[10] + x[6]*y[8] + x[8]*y[6] - x[10]*y[4] + x[12]*y[2] + x[14]*y[0] - x[16]*y[30] + x[18]*y[28] - x[20]*y[26] - x[22]*y[24] + x[24]*y[22] + x[26]*y[20] - x[28]*y[18] + x[30]*y[16],
        x[0]*y[15] + x[1]*y[14] - x[2]*y[13] + x[3]*y[12] + x[4]*y[11] - x[5]*y[10] + x[6]*y[9] + x[7]*y[8] - x[8]*y[7] + x[9]*y[6] - x[10]*y[5] - x[11]*y[4] + x[12]*y[3] + x[13]*y[2] - x[14]*y[1] + x[15]*y[0] + x[16]*y[31] - x[17]*y[30] + x[18]*y[29] + x[19]*y[28] - x[20]*y[27] - x[21]*y[26] + x[22]*y[25] - x[23]*y[24] + x[24]*y[23] + x[25]*y[22] - x[26]*y[21] + x[27]*y[20] + x[28]*y[19] - x[29]*y[18] + x[30]*y[17] + x[31]*y[16],
        x[0]*y[16] + x[2]*y[18] + x[4]*y[20] - x[6]*y[22] + x[8]*y[24] - x[10]*y[26] - x[12]*y[28] - x[14]*y[30] + x[16]*y[0] - x[18]*y[2] - x[20]*y[4] - x[22]*y[6] - x[24]*y[8] - x[26]*y[10] - x[28]*y[12] + x[30]*y[14],
        x[0]*y[17] + x[1]*y[16] - x[2]*y[19] + x[3]*y[18] - x[4]*y[21] + x[5]*y[20] - x[6]*y[23] - x[7]*y[22] - x[8]*y[25] + x[9]*y[24] - x[10]*y[27] - x[11]*y[26] - x[12]*y[29] - x[13]*y[28] + x[14]*y[31] - x[15]*y[30] - x[16]*y[1] + x[17]*y[0] - x[18]*y[3] - x[19]*y[2] - x[20]*y[5] - x[21]*y[4] + x[22]*y[7] - x[23]*y[6] - x[24]*y[9] - x[25]*y[8] + x[26]*y[11] - x[27]*y[10] + x[28]*y[13] - x[29]*y[12] + x[30]*y[15] + x[31]*y[14],
        x[0]*y[18] + x[2]*y[16] - x[4]*y[22] + x[6]*y[20] - x[8]*y[26] + x[10]*y[24] - x[12]*y[30] - x[14]*y[28] - x[16]*y[2] + x[18]*y[0] - x[20]*y[6] - x[22]*y[4] - x[24]*y[10] - x[26]*y[8] + x[28]*y[14] - x[30]*y[12],
        x[0]*y[19] + x[1]*y[18] - x[2]*y[17] + x[3]*y[16] + x[4]*y[23] - x[5]*y[22] + x[6]*y[21] + x[7]*y[20] + x[8]*y[27] - x[9]*y[26] + x[10]*y[25] + x[11]*y[24] - x[12]*y[31] - x[13]*y[30] + x[14]*y[29] - x[15]*y[28] + x[16]*y[3] - x[17]*y[2] + x[18]*y[1] + x[19]*y[0] - x[20]*y[7] - x[21]*y[6] + x[22]*y[5] - x[23]*y[4] - x[24]*y[11] - x[25]*y[10] + x[26]*y[9] - x[27]*y[8] - x[28]*y[15] + x[29]*y[14] - x[30]*y[13] - x[31]*y[12],
        x[0]*y[20] + x[2]*y[22] + x[4]*y[16] - x[6]*y[18] - x[8]*y[28] + x[10]*y[30] + x[12]*y[24] + x[14]*y[26] - x[16]*y[4] + x[18]*y[6] + x[20]*y[0] + x[22]*y[2] - x[24]*y[12] - x[26]*y[14] - x[28]*y[8] + x[30]*y[10],
        x[0]*y[21] + x[1]*y[20] - x[2]*y[23] + x[3]*y[22] - x[4]*y[17] + x[5]*y[16] - x[6]*y[19] - x[7]*y[18] + x[8]*y[29] - x[9]*y[28] + x[10]*y[31] + x[11]*y[30] + x[12]*y[25] + x[13]*y[24] - x[14]*y[27] + x[15]*y[26] + x[16]*y[5] - x[17]*y[4] + x[18]*y[7] + x[19]*y[6] + x[20]*y[1] + x[21]*y[0] - x[22]*y[3] + x[23]*y[2] - x[24]*y[13] - x[25]*y[12] + x[26]*y[15] - x[27]*y[14] + x[28]*y[9] - x[29]*y[8] + x[30]*y[11] + x[31]*y[10],
        x[0]*y[22] + x[2]*y[20] - x[4]*y[18] + x[6]*y[16] + x[8]*y[30] - x[10]*y[28] + x[12]*y[26] + x[14]*y[24] + x[16]*y[6] - x[18]*y[4] + x[20]*y[2] + x[22]*y[0] - x[24]*y[14] - x[26]*y[12] + x[28]*y[10] - x[30]*y[8],
        x[0]*y[23] + x[1]*y[22] - x[2]*y[21] + x[3]*y[20] + x[4]*y[19] - x[5]*y[18] + x[6]*y[17] + x[7]*y[16] - x[8]*y[31] + x[9]*y[30] - x[10]*y[29] - x[11]*y[28] + x[12]*y[27] + x[13]*y[26] - x[14]*y[25] + x[15]*y[24] - x[16]*y[7] + x[17]*y[6] - x[18]*y[5] - x[19]*y[4] + x[20]*y[3] + x[21]*y[2] - x[22]*y[1] + x[23]*y[0] - x[24]*y[15] - x[25]*y[14] + x[26]*y[13] - x[27]*y[12] - x[28]*y[11] + x[29]*y[10] - x[30]*y[9] - x[31]*y[8],
        x[0]*y[24] + x[2]*y[26] + x[4]*y[28] - x[6]*y[30] + x[8]*y[16] - x[10]*y[18] - x[12]*y[20] - x[14]*y[22] - x[16]*y[8] + x[18]*y[10] + x[20]*y[12] + x[22]*y[14] + x[24]*y[0] + x[26]*y[2] + x[28]*y[4] - x[30]*y[6],
        x[0]*y[25] + x[1]*y[24] - x[2]*y[27] + x[3]*y[26] - x[4]*y[29] + x[5]*y[28] - x[6]*y[31] - x[7]*y[30] - x[8]*y[17] + x[9]*y[16] - x[10]*y[19] - x[11]*y[18] - x[12]*y[21] - x[13]*y[20] + x[14]*y[23] - x[15]*y[22] + x[16]*y[9] - x[17]*y[8] + x[18]*y[11] + x[19]*y[10] + x[20]*y[13] + x[21]*y[12] - x[22]*y[15] + x[23]*y[14] + x[24]*y[1] + x[25]*y[0] - x[26]*y[3] + x[27]*y[2] - x[28]*y[5] + x[29]*y[4] - x[30]*y[7] - x[31]*y[6],
        x[0]*y[26] + x[2]*y[24] - x[4]*y[30] + x[6]*y[28] - x[8]*y[18] + x[10]*y[16] - x[12]*y[22] - x[14]*y[20] + x[16]*y[10] - x[18]*y[8] + x[20]*y[14] + x[22]*y[12] + x[24]*y[2] + x[26]*y[0] - x[28]*y[6] + x[30]*y[4],
        x[0]*y[27] + x[1]*y[26] - x[2]*y[25] + x[3]*y[24] + x[4]*y[31] - x[5]*y[30] + x[6]*y[29] + x[7]*y[28] + x[8]*y[19] - x[9]*y[18] + x[10]*y[17] + x[11]*y[16] - x[12]*y[23] - x[13]*y[22] + x[14]*y[21] - x[15]*y[20] - x[16]*y[11] + x[17]*y[10] - x[18]*y[9] - x[19]*y[8] + x[20]*y[15] + x[21]*y[14] - x[22]*y[13] + x[23]*y[12] + x[24]*y[3] + x[25]*y[2] - x[26]*y[1] + x[27]*y[0] + x[28]*y[7] - x[29]*y[6] + x[30]*y[5] + x[31]*y[4],
        x[0]*y[28] + x[2]*y[30] + x[4]*y[24] - x[6]*y[26] - x[8]*y[20] + x[10]*y[22] + x[12]*y[16] + x[14]*y[18] + x[16]*y[12] - x[18]*y[14] - x[20]*y[8] - x[22]*y[10] + x[24]*y[4] + x[26]*y[6] + x[28]*y[0] - x[30]*y[2],
        x[0]*y[29] + x[1]*y[28] - x[2]*y[31] + x[3]*y[30] - x[4]*y[25] + x[5]*y[24] - x[6]*y[27] - x[7]*y[26] + x[8]*y[21] - x[9]*y[20] + x[10]*y[23] + x[11]*y[22] + x[12]*y[17] + x[13]*y[16] - x[14]*y[19] + x[15]*y[18] - x[16]*y[13] + x[17]*y[12] - x[18]*y[15] - x[19]*y[14] - x[20]*y[9] - x[21]*y[8] + x[22]*y[11] - x[23]*y[10] + x[24]*y[5] + x[25]*y[4] - x[26]*y[7] + x[27]*y[6] - x[28]*y[1] + x[29]*y[0] - x[30]*y[3] - x[31]*y[2],
        x[0]*y[30] + x[2]*y[28] - x[4]*y[26] + x[6]*y[24] + x[8]*y[22] - x[10]*y[20] + x[12]*y[18] + x[14]*y[16] - x[16]*y[14] + x[18]*y[12] - x[20]*y[10] - x[22]*y[8] + x[24]*y[6] + x[26]*y[4] - x[28]*y[2] + x[30]*y[0],
        x[0]*y[31] + x[1]*y[30] - x[2]*y[29] + x[3]*y[28] + x[4]*y[27] - x[5]*y[26] + x[6]*y[25] + x[7]*y[24] - x[8]*y[23] + x[9]*y[22] - x[10]*y[21] - x[11]*y[20] + x[12]*y[19] + x[13]*y[18] - x[14]*y[17] + x[15]*y[16] + x[16]*y[15] - x[17]*y[14] + x[18]*y[13] + x[19]*y[12] - x[20]*y[11] - x[21]*y[10] + x[22]*y[9] - x[23]*y[8] + x[24]*y[7] + x[25]*y[6] - x[26]*y[5] + x[27]*y[4] + x[28]*y[3] - x[29]*y[2] + x[30]*y[1] + x[31]*y[0]});
}

static inline pga4
pga4_square(pga4 a)
{
    const double *x = a.data;
    return pga4_make((const double[]){
        x[0]*x[0] + x[2]*x[2] + x[4]*x[4] - x[6]*x[6] + x[8]*x[8] - x[10]*x[10] - x[12]*x[12] - x[14]*x[14] + x[16]*x[16] - x[18]*x[18] - x[20]*x[20] - x[22]*x[22] - x[24]*x[24] - x[26]*x[26] - x[28]*x[28] + x[30]*x[30],
        2*x[0]*x[1] - 2*x[6]*x[7] - 2*x[10]*x[11] - 2*x[12]*x[13] - 2*x[18]*x[19] - 2*x[20]*x[21] - 2*x[24]*x[25] + 2*x[30]*x[31],
        2*x[0]*x[2] - 2*x[12]*x[14] - 2*x[20]*x[22] - 2*x[24]*x[26],
        2*x[0]*x[3] + 2*x[4]*x[7] + 2*x[8]*x[11] - 2*x[12]*x[15] + 2*x[16]*x[19] - 2*x[20]*x[23] - 2*x[24]*x[27] - 2*x[28]*x[31],
        2*x[0]*x[4] + 2*x[10]*x[14] + 2*x[18]*x[22] - 2*x[24]*x[28],
        2*x[0]*x[5] - 2*x[2]*x[7] + 2*x[8]*x[13] + 2*x[10]*x[15] + 2*x[16]*x[21] + 2*x[18]*x[23] - 2*x[24]*x[29] + 2*x[26]*x[31],
        2*x[0]*x[6] + 2*x[8]*x[14] + 2*x[16]*x[22] - 2*x[24]*x[30],
        2*x[0]*x[7] + 2*x[1]*x[6] - 2*x[2]*x[5] + 2*x[3]*x[4] - 2*x[24]*x[31] - 2*x[25]*x[30] + 2*x[26]*x[29] - 2*x[27]*x[28],
        2*x[0]*x[8] - 2*x[6]*x[14] + 2*x[18]*x[26] + 2*x[20]*x[28],
        2*x[0]*x[9] - 2*x[2]*x[11] - 2*x[4]*x[13] - 2*x[6]*x[15] + 2*x[16]*x[25] + 2*x[18]*x[27] + 2*x[20]*x[29] - 2*x[22]*x[31],
        2*x[0]*x[10] - 2*x[4]*x[14] + 2*x[16]*x[26] + 2*x[20]*x[30],
        2*x[0]*x[11] + 2*x[1]*x[10] - 2*x[2]*x[9] + 2*x[3]*x[8] + 2*x[20]*x[31] + 2*x[21]*x[30] - 2*x[22]*x[29] + 2*x[23]*x[28],
        2*x[0]*x[12] + 2*x[2]*x[14] + 2*x[16]*x[28] - 2*x[18]*x[30],
        2*x[0]*x[13] + 2*x[1]*x[12] - 2*x[4]*x[9] + 2*x[5]*x[8] - 2*x[18]*x[31] - 2*x[19]*x[30] + 2*x[22]*x[27] - 2*x[23]*x[26],
        2*x[0]*x[14] + 2*x[2]*x[12] - 2*x[4]*x[10] + 2*x[6]*x[8],
        2*x[0]*x[15] + 2*x[3]*x[12] - 2*x[5]*x[10] + 2*x[6]*x[9] + 2*x[16]*x[31] + 2*x[19]*x[28] - 2*x[21]*x[26] + 2*x[22]*x[25],
        2*x[0]*x[16] - 2*x[6]*x[22] - 2*x[10]*x[26] - 2*x[12]*x[28],
        2*x[0]*x[17] - 2*x[2]*x[19] - 2*x[4]*x[21] - 2*x[6]*x[23] - 2*x[8]*x[25] - 2*x[10]*x[27] - 2*x[12]*x[29] + 2*x[14]*x[31],
        2*x[0]*x[18] - 2*x[4]*x[22] - 2*x[8]*x[26] - 2*x[12]*x[30],
        2*x[0]*x[19] + 2*x[1]*x[18] - 2*x[2]*x[17] + 2*x[3]*x[16] - 2*x[12]*x[31] - 2*x[13]*x[30] + 2*x[14]*x[29] - 2*x[15]*x[28],
        2*x[0]*x[20] + 2*x[2]*x[22] - 2*x[8]*x[28] + 2*x[10]*x[30],
        2*x[0]*x[21] + 2*x[1]*x[20] - 2*x[4]*x[17] + 2*x[5]*x[16] + 2*x[10]*x[31] + 2*x[11]*x[30] - 2*x[14]*x[27] + 2*x[15]*x[26],
        2*x[0]*x[22] + 2*x[2]*x[20] - 2*x[4]*x[18] + 2*x[6]*x[16],
        2*x[0]*x[23] + 2*x[3]*x[20] - 2*x[5]*x[18] + 2*x[6]*x[17] - 2*x[8]*x[31] - 2*x[11]*x[28] + 2*x[13]*x[26] - 2*x[14]*x[25],
        2*x[0]*x[24] + 2*x[2]*x[26] + 2*x[4]*x[28] - 2*x[6]*x[30],
        2*x[0]*x[25] + 2*x[1]*x[24] - 2*x[6]*x[31] - 2*x[7]*x[30] - 2*x[8]*x[17] + 2*x[9]*x[16] + 2*x[14]*x[23] - 2*x[15]*x[22],
        2*x[0]*x[26] + 2*x[2]*x[24] - 2*x[8]*x[18] + 2*x[10]*x[16],
        2*x[0]*x[27] + 2*x[3]*x[24] + 2*x[4]*x[31] + 2*x[7]*x[28] - 2*x[9]*x[18] + 2*x[10]*x[17] - 2*x[13]*x[22] + 2*x[14]*x[21],
        2*x[0]*x[28] + 2*x[4]*x[24] - 2*x[8]*x[20] + 2*x[12]*x[16],
        2*x[0]*x[29] - 2*x[2]*x[31] + 2*x[5]*x[24] - 2*x[7]*x[26] - 2*x[9]*x[20] + 2*x[11]*x[22] + 2*x[12]*x[17] - 2*x[14]*x[19],
        2*x[0]*x[30] + 2*x[6]*x[24] - 2*x[10]*x[20] + 2*x[12]*x[18],
        2*x[0]*x[31] + 2*x[1]*x[30] - 2*x[2]*x[29] + 2*x[3]*x[28] + 2*x[4]*x[27] - 2*x[5]*x[26] + 2*x[6]*x[25] + 2*x[7]*x[24] - 2*x[8]*x[23] + 2*x[9]*x[22] - 2*x[10]*x[21] - 2*x[11]*x[20] + 2*x[12]*x[19] + 2*x[13]*x[18] - 2*x[14]*x[17] + 2*x[15]*x[16]});
}

static inline pga4
pga4_conjugate(pga4 a)
{
    return pga4_make((const double[]){a.data[0], -a.data[1], -a.data[2], -a.data[3], -a.data[4], -a.data[5], -a.data[6], a.data[7], -a.data[8], -a.data[9], -a.data[10], a.data[11], -a.data[12], a.data[13], a.data[14], a.data[15], -a.data[16], -a.data[17], -a.data[18], a.data[19], -a.data[20], a.data[21], a.data[22], a.data[23], -a.data[24], a.data[25], a.data[26], a.data[27], a.data[28], a.data[29], a.data[30], -a.data[31]});
}

static inline void
pga4_basis(pga4 result[32])
{
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            result[i].data[j] = i == j ? 1 : 0;
        }
    }
}

static inline void
pga4_gens(pga4 result[5])
{
    result[0] = pga4_make((const double[ALGPY_PGA4_N]){[1] = 1});
    result[1] = pga4_make((const double[ALGPY_PGA4_N]){[2] = 1});
    result[2] = pga4_make((const double[ALGPY_PGA4_N]){[4] = 1});
    result[3] = pga4_make((const double[ALGPY_PGA4_N]){[8] = 1});
    result[4] = pga4_make((const double[ALGPY_PGA4_N]){[16] = 1});
}

static inline void
pga4_coeffs(pga4 result[1024])
{
    pga4 basis[32];
    int i, j;
    pga4_basis(basis);
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            result[i*32+j] = pga4_mul(basis[i], basis[j]);
        }
    }
}

static inline double
pga4_scalar_product(pga4 a, pga4 b)
{
    const double *x = a.data, *y = b.data;
    return x[0]*y[0] + x[2]*y[2] + x[4]*y[4] - x[6]*y[6] + x[8]*y[8] - x[10]*y[10] - x[12]*y[12] - x[14]*y[14] + x[16]*y[16] - x[18]*y[18] - x[20]*y[20] - x[22]*y[22] - x[24]*y[24] - x[26]*y[26] - x[28]*y[28] + x[30]*y[30];
}

/* The Euclidean quotient Cl(4,0) has a real quartic. The square-zero
 * radical can double its multiplicities, giving a generic degree of eight. */
static inline void
pga4_quotient_minpoly(pga4 a, double result[5])
{
    double w = a.data[0], s2, s3, s4;
    pga4 b = a, b2;
    b.data[0] = 0;
    b2 = pga4_square(b);
    s2 = b2.data[0];
    s3 = pga4_scalar_product(b2, b);
    s4 = pga4_scalar_product(b2, b2);
    result[0] = w*w*w*w - 2*w*w*s2 + 4*w*s3/3 + 2*s2*s2 - s4;
    result[1] = -4*w*w*w + 4*w*s2 - 4*s3/3;
    result[2] = 6*w*w - 2*s2;
    result[3] = -4*w;
    result[4] = 1;
}

static inline void
pga4_minpoly(pga4 a, double result[9])
{
    double p[5];
    pga4_quotient_minpoly(a, p);
    for (int i = 0; i < 9; i++) result[i] = 0;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) result[i+j] += p[i]*p[j];
    }
}

static inline double
pga4_norm8(pga4 a)
{
    double p[5];
    pga4_quotient_minpoly(a, p);
    return p[0]*p[0];
}

static inline double
pga4_norm(pga4 a)
{
    double p[5];
    pga4_quotient_minpoly(a, p);
    return sqrt(sqrt(fabs(p[0])));
}

static inline void
pga4_eivals(pga4 a, complex result[8])
{
    double w = a.data[0], p[5];
    complex q[5];
    a.data[0] = 0;
    pga4_quotient_minpoly(a, p);
    for (int i = 0; i < 5; i++) q[i] = complex_from_re(p[i]);
    poly4_complex_coeffs_to_roots(q, result);
    for (int i = 0; i < 4; i++) {
        result[i].data[0] += w;
        result[i+4] = result[i];
    }
}

static inline pga4
pga4_inv(pga4 a)
{
    double p[5];
    pga4 u = a;
    for (int i = 1; i < 32; i += 2) u.data[i] = 0;
    pga4_quotient_minpoly(u, p);
    pga4 u2 = pga4_square(u), u3 = pga4_mul(u2, u);
    pga4 inverse = vec32_div_scalar(vec32_neg(vec32_add(u3,
        vec32_add(vec32_scale(u2, p[3]), vec32_add(vec32_scale(u, p[2]), pga4_from_re(p[1]))))), p[0]);
    pga4 radical = vec32_sub(a, u);
    return vec32_sub(inverse, pga4_mul(pga4_mul(inverse, radical), inverse));
}

static inline pga4
pga4_div(pga4 a, pga4 b)
{
    return pga4_mul(a, pga4_inv(b));
}

static inline pga4
pga4_nan(void)
{
    pga4 result;
    for (int i = 0; i < 32; i++) result.data[i] = NAN;
    return result;
}

static inline pga4
pga4_real_slice(complex value)
{
    pga4 result = pga4_from_re(value.data[0]);
    result.data[6] = value.data[1];
    return result;
}

static inline int
pga4_unused_generator(pga4 a)
{
    int generator, i;
    for (generator = 1; generator < 5; generator++) {
        int used = 0;
        for (i = 1; i < 32; i++) {
            if ((i & (1 << generator)) && a.data[i] != 0) {
                used = 1;
                break;
            }
        }
        if (!used) {
            return generator;
        }
    }
    return -1;
}

static inline pga4
pga4_apply_subalgebra(pga4 a, int generator,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    double (*deriv4)(double), double (*deriv5)(double),
    double (*deriv6)(double), double (*deriv7)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex),
    complex (*cderiv2)(complex), complex (*cderiv3)(complex))
{
    (void)deriv4; (void)deriv5; (void)deriv6; (void)deriv7;
    (void)cderiv2; (void)cderiv3;
    pga3 suba;
    pga4 result = {{0}};
    int low = (1 << generator)-1;
    for (int i = 0; i < 16; i++) suba.data[i] = a.data[(i & low) | ((i & ~low) << 1)];
    pga3 value = pga3_apply(suba, f, deriv1, deriv2, deriv3, cf, cderiv1);
    if (isnan(value.data[0])) return pga4_nan();
    for (int i = 0; i < 16; i++) result.data[(i & low) | ((i & ~low) << 1)] = value.data[i];
    return result;
}

static inline int
pga4_quadratic_relation(pga4 b, pga4 b2, double *linear)
{
    int i, pivot = 1;
    double bound, size = 0;
    for (i = 2; i < 32; i++) {
        if (fabs(b.data[i]) > fabs(b.data[pivot])) {
            pivot = i;
        }
    }
    *linear = b2.data[pivot]/b.data[pivot];
    if (!isfinite(*linear)) {
        return 0;
    }
    for (i = 0; i < 32; i++) {
        size = fmax(size, fabs(b2.data[i]));
    }
    bound = 32*DBL_EPSILON*(size + fabs(*linear*b.data[pivot]));
    for (i = 1; i < 32; i++) {
        double value = *linear*b.data[i];
        if (fabs(b2.data[i] - value) > bound) {
            return 0;
        }
    }
    return 1;
}

static inline int
pga4_cubic_relation(pga4 b, pga4 b2, pga4 b3, double p[4])
{
    int i, first = 1, second = 1;
    double ratio, largest = 0, linear, quadratic, bound = 0;
    for (i = 2; i < 32; i++) {
        if (fabs(b.data[i]) > fabs(b.data[first])) {
            first = i;
        }
    }
    ratio = b2.data[first]/b.data[first];
    for (i = 1; i < 32; i++) {
        double residual = fabs(b2.data[i] - ratio*b.data[i]);
        if (residual > largest) {
            largest = residual;
            second = i;
        }
    }
    if (largest == 0) {
        return 0;
    }
    quadratic = (b3.data[second] - b3.data[first]/b.data[first]*b.data[second])/
        (b2.data[second] - ratio*b.data[second]);
    linear = (b3.data[first] - quadratic*b2.data[first])/b.data[first];
    if (!isfinite(linear) || !isfinite(quadratic)) {
        return 0;
    }
    for (i = 0; i < 32; i++) {
        bound = fmax(bound, fabs(b3.data[i]) + fabs(quadratic*b2.data[i]) + fabs(linear*b.data[i]));
    }
    bound *= 64*DBL_EPSILON;
    for (i = 1; i < 32; i++) {
        double q = quadratic*b2.data[i], l = linear*b.data[i];
        if (fabs(b3.data[i] - q - l) > bound) {
            return 0;
        }
    }
    p[0] = quadratic*b2.data[0] - b3.data[0];
    p[1] = -linear;
    p[2] = -quadratic;
    p[3] = 1;
    return 1;
}

/* Check whether a quotient quartic also annihilates the radical part.
 * This avoids derivatives at semisimple zero eigenvalues. */
static inline int
pga4_quartic_relation(pga4 b, pga4 b2, pga4 b3, const double p[5])
{
    pga4 b4 = pga4_square(b2);
    double error = 0, bound = fabs(p[0]);
    for (int i = 0; i < 32; i++) {
        double value = b4.data[i] + p[3]*b3.data[i] + p[2]*b2.data[i] + p[1]*b.data[i] + (i == 0 ? p[0] : 0);
        error = fmax(error, fabs(value));
        bound = fmax(bound, fabs(b4.data[i]) + fabs(p[3]*b3.data[i]) + fabs(p[2]*b2.data[i]) + fabs(p[1]*b.data[i]));
    }
    return error <= 64*DBL_EPSILON*bound;
}

static inline pga4
pga4_apply(pga4 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    double (*deriv4)(double), double (*deriv5)(double),
    double (*deriv6)(double), double (*deriv7)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex),
    complex (*cderiv2)(complex), complex (*cderiv3)(complex))
{
    double w = a.data[0], linear, coefficients[8] = {0}, cubic[4];
    if (vec32_is_real(a)) return pga4_real_slice(cf(complex_from_re(w)));
    int generator = pga4_unused_generator(a);
    if (generator >= 0) return pga4_apply_subalgebra(a, generator,
        f, deriv1, deriv2, deriv3, deriv4, deriv5, deriv6, deriv7, cf, cderiv1, cderiv2, cderiv3);
    pga4 b = a;
    b.data[0] = 0;
    pga4 b2 = pga4_square(b), b3;
    int degree;
    if (pga4_quadratic_relation(b, b2, &linear)) {
        double shift = linear/2;
        complex root = complex_sqrt_real(b2.data[0] + shift*shift);
        complex nodes[2] = {complex_ladd(vec2_neg(root), shift), complex_ladd(root, shift)};
        if (!hermite_interpolate(2, w, nodes, f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients)) return pga4_nan();
        pga4 result = vec32_scale(b, coefficients[1]);
        result.data[0] += coefficients[0];
        return result;
    }
    b3 = pga4_mul(b2, b);
    if (pga4_cubic_relation(b, b2, b3, cubic)) {
        complex nodes[3];
        poly3_coeffs_to_roots(cubic, nodes);
        if (!hermite_interpolate(3, w, nodes, f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients)) return pga4_nan();
        degree = 3;
    } else {
        double s2 = b2.data[0], s3 = b3.data[0], s4 = pga4_scalar_product(b2, b2);
        double q[5] = {2*s2*s2 - s4, -4*s3/3, -2*s2, 0, 1};
        if (pga4_quartic_relation(b, b2, b3, q)) {
            if (!poly4_apply(q, w, f, deriv1, deriv2, deriv3, cf, cderiv1, coefficients)) return pga4_nan();
            degree = 4;
        } else {
            complex nodes[8], p[5];
            for (int i = 0; i < 5; i++) p[i] = complex_from_re(q[i]);
            poly4_complex_coeffs_to_roots(p, nodes);
            for (int i = 0; i < 4; i++) nodes[4+i] = nodes[i];
            double (*derivatives[7])(double) = {deriv1, deriv2, deriv3, deriv4, deriv5, deriv6, deriv7};
            complex (*cderivatives[3])(complex) = {cderiv1, cderiv2, cderiv3};
            if (!hermite8_interpolate(8, w, nodes, f, derivatives, cf, cderivatives, coefficients)) return pga4_nan();
            degree = 8;
        }
    }
    pga4 result, power = b3;
    for (int i = 0; i < 32; i++) result.data[i] = coefficients[1]*b.data[i] + coefficients[2]*b2.data[i] + coefficients[3]*b3.data[i];
    result.data[0] += coefficients[0];
    for (int i = 4; i < degree; i++) {
        power = pga4_mul(power, b);
        result = vec32_add(result, vec32_scale(power, coefficients[i]));
    }
    return result;
}

static inline pga4
pga4_exp(pga4 a)
{
    return pga4_apply(a, double_exp,
        double_exp_deriv1, double_exp_deriv2, double_exp_deriv3, double_exp_deriv4, double_exp_deriv5, double_exp_deriv6, double_exp_deriv7,
        complex_exp, complex_exp_deriv1, complex_exp_deriv2, complex_exp_deriv3);
}

static inline pga4
pga4_expm1(pga4 a)
{
    return pga4_apply(a, double_expm1,
        double_expm1_deriv1, double_expm1_deriv2, double_expm1_deriv3, double_expm1_deriv4, double_expm1_deriv5, double_expm1_deriv6, double_expm1_deriv7,
        complex_expm1, complex_expm1_deriv1, complex_expm1_deriv2, complex_expm1_deriv3);
}

static inline pga4
pga4_log(pga4 a)
{
    return pga4_apply(a, double_log,
        double_log_deriv1, double_log_deriv2, double_log_deriv3, double_log_deriv4, double_log_deriv5, double_log_deriv6, double_log_deriv7,
        complex_log, complex_log_deriv1, complex_log_deriv2, complex_log_deriv3);
}

static inline pga4
pga4_log1p(pga4 a)
{
    return pga4_apply(a, double_log1p,
        double_log1p_deriv1, double_log1p_deriv2, double_log1p_deriv3, double_log1p_deriv4, double_log1p_deriv5, double_log1p_deriv6, double_log1p_deriv7,
        complex_log1p, complex_log1p_deriv1, complex_log1p_deriv2, complex_log1p_deriv3);
}

static inline pga4
pga4_sqrt(pga4 a)
{
    return pga4_apply(a, double_sqrt,
        double_sqrt_deriv1, double_sqrt_deriv2, double_sqrt_deriv3, double_sqrt_deriv4, double_sqrt_deriv5, double_sqrt_deriv6, double_sqrt_deriv7,
        complex_sqrt, complex_sqrt_deriv1, complex_sqrt_deriv2, complex_sqrt_deriv3);
}

static inline pga4
pga4_sin(pga4 a)
{
    return pga4_apply(a, double_sin,
        double_sin_deriv1, double_sin_deriv2, double_sin_deriv3, double_sin_deriv4, double_sin_deriv5, double_sin_deriv6, double_sin_deriv7,
        complex_sin, complex_sin_deriv1, complex_sin_deriv2, complex_sin_deriv3);
}

static inline pga4
pga4_cos(pga4 a)
{
    return pga4_apply(a, double_cos,
        double_cos_deriv1, double_cos_deriv2, double_cos_deriv3, double_cos_deriv4, double_cos_deriv5, double_cos_deriv6, double_cos_deriv7,
        complex_cos, complex_cos_deriv1, complex_cos_deriv2, complex_cos_deriv3);
}

static inline pga4
pga4_tan(pga4 a)
{
    return pga4_apply(a, double_tan,
        double_tan_deriv1, double_tan_deriv2, double_tan_deriv3, double_tan_deriv4, double_tan_deriv5, double_tan_deriv6, double_tan_deriv7,
        complex_tan, complex_tan_deriv1, complex_tan_deriv2, complex_tan_deriv3);
}

static inline pga4
pga4_asin(pga4 a)
{
    return pga4_apply(a, double_asin,
        double_asin_deriv1, double_asin_deriv2, double_asin_deriv3, double_asin_deriv4, double_asin_deriv5, double_asin_deriv6, double_asin_deriv7,
        complex_asin, complex_asin_deriv1, complex_asin_deriv2, complex_asin_deriv3);
}

static inline pga4
pga4_acos(pga4 a)
{
    return pga4_apply(a, double_acos,
        double_acos_deriv1, double_acos_deriv2, double_acos_deriv3, double_acos_deriv4, double_acos_deriv5, double_acos_deriv6, double_acos_deriv7,
        complex_acos, complex_acos_deriv1, complex_acos_deriv2, complex_acos_deriv3);
}

static inline pga4
pga4_atan(pga4 a)
{
    return pga4_apply(a, double_atan,
        double_atan_deriv1, double_atan_deriv2, double_atan_deriv3, double_atan_deriv4, double_atan_deriv5, double_atan_deriv6, double_atan_deriv7,
        complex_atan, complex_atan_deriv1, complex_atan_deriv2, complex_atan_deriv3);
}

static inline pga4
pga4_sinh(pga4 a)
{
    return pga4_apply(a, double_sinh,
        double_sinh_deriv1, double_sinh_deriv2, double_sinh_deriv3, double_sinh_deriv4, double_sinh_deriv5, double_sinh_deriv6, double_sinh_deriv7,
        complex_sinh, complex_sinh_deriv1, complex_sinh_deriv2, complex_sinh_deriv3);
}

static inline pga4
pga4_cosh(pga4 a)
{
    return pga4_apply(a, double_cosh,
        double_cosh_deriv1, double_cosh_deriv2, double_cosh_deriv3, double_cosh_deriv4, double_cosh_deriv5, double_cosh_deriv6, double_cosh_deriv7,
        complex_cosh, complex_cosh_deriv1, complex_cosh_deriv2, complex_cosh_deriv3);
}

static inline pga4
pga4_tanh(pga4 a)
{
    return pga4_apply(a, double_tanh,
        double_tanh_deriv1, double_tanh_deriv2, double_tanh_deriv3, double_tanh_deriv4, double_tanh_deriv5, double_tanh_deriv6, double_tanh_deriv7,
        complex_tanh, complex_tanh_deriv1, complex_tanh_deriv2, complex_tanh_deriv3);
}

static inline pga4
pga4_asinh(pga4 a)
{
    return pga4_apply(a, double_asinh,
        double_asinh_deriv1, double_asinh_deriv2, double_asinh_deriv3, double_asinh_deriv4, double_asinh_deriv5, double_asinh_deriv6, double_asinh_deriv7,
        complex_asinh, complex_asinh_deriv1, complex_asinh_deriv2, complex_asinh_deriv3);
}

static inline pga4
pga4_acosh(pga4 a)
{
    return pga4_apply(a, double_acosh,
        double_acosh_deriv1, double_acosh_deriv2, double_acosh_deriv3, double_acosh_deriv4, double_acosh_deriv5, double_acosh_deriv6, double_acosh_deriv7,
        complex_acosh, complex_acosh_deriv1, complex_acosh_deriv2, complex_acosh_deriv3);
}

static inline pga4
pga4_atanh(pga4 a)
{
    return pga4_apply(a, double_atanh,
        double_atanh_deriv1, double_atanh_deriv2, double_atanh_deriv3, double_atanh_deriv4, double_atanh_deriv5, double_atanh_deriv6, double_atanh_deriv7,
        complex_atanh, complex_atanh_deriv1, complex_atanh_deriv2, complex_atanh_deriv3);
}

static inline pga4
pga4_radd(double a, pga4 b)
{
    return pga4_add(pga4_from_re(a), b);
}

static inline pga4
pga4_ladd(pga4 a, double b)
{
    return pga4_add(a, pga4_from_re(b));
}

static inline pga4
pga4_rsub(double a, pga4 b)
{
    return pga4_sub(pga4_from_re(a), b);
}

static inline pga4
pga4_lsub(pga4 a, double b)
{
    return pga4_sub(a, pga4_from_re(b));
}

static inline pga4
pga4_rmul(double a, pga4 b)
{
    return pga4_mul(pga4_from_re(a), b);
}

static inline pga4
pga4_lmul(pga4 a, double b)
{
    return pga4_mul(a, pga4_from_re(b));
}

static inline pga4
pga4_rdiv(double a, pga4 b)
{
    return pga4_div(pga4_from_re(a), b);
}

static inline pga4
pga4_ldiv(pga4 a, double b)
{
    return pga4_div(a, pga4_from_re(b));
}

static inline pga4
pga4_ipow(pga4 a, int64_t b)
{
    pga4 result = pga4_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return pga4_mul(a, a);
    }
    if (b == -1) {
        return pga4_inv(a);
    }
    if (b < 0) {
        a = pga4_inv(a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = pga4_mul(result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = pga4_mul(a, a);
        }
    }
    return result;
}

static inline pga4
pga4_rpow(pga4 a, pga4 b)
{
    return pga4_exp(pga4_mul(pga4_log(a), b));
}

static inline pga4
pga4_pow(pga4 a, pga4 b)
{
    return pga4_rpow(a, b);
}

static inline pga4
pga4_lpow(pga4 a, pga4 b)
{
    return pga4_exp(pga4_mul(b, pga4_log(a)));
}

static inline pga4
pga4_rpows(double a, pga4 b)
{
    return pga4_rpow(pga4_from_re(a), b);
}

static inline pga4
pga4_lpows(pga4 a, double b)
{
    return pga4_rpow(a, pga4_from_re(b));
}

static inline pga4
pga4_mix(pga4 a, pga4 b, double t)
{
    pga4 delta = pga4_mul(pga4_inv(a), b);
    return pga4_mul(a, pga4_lpows(delta, t));
}

static inline void
pga4_lmat(pga4 a, double result[ALGPY_PGA4_N*ALGPY_PGA4_N])
{
    pga4 basis = {{0}};
    pga4 product;
    int row, col;
    for (col = 0; col < ALGPY_PGA4_N; col++) {
        basis.data[col] = 1;
        product = pga4_mul(a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_PGA4_N; row++) {
            result[row * ALGPY_PGA4_N + col] = product.data[row];
        }
    }
}

static inline void
pga4_rmat(pga4 a, double result[ALGPY_PGA4_N*ALGPY_PGA4_N])
{
    pga4 basis = {{0}};
    pga4 product;
    int row, col;
    for (row = 0; row < ALGPY_PGA4_N; row++) {
        basis.data[row] = 1;
        product = pga4_mul(basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_PGA4_N; col++) {
            result[row * ALGPY_PGA4_N + col] = product.data[col];
        }
    }
}

/* Lattice API. The exterior and regressive products and the Poincare duals
 * do not depend on the metric, so they are those of exterior5, which shares
 * this coefficient order: join with identity bottom = 1, meet with identity
 * top = the pseudoscalar, rdual(e_a) join e_a = top = e_a join ldual(e_a). */
static inline pga4
pga4_bottom(void)
{
    return exterior5_bottom();
}

static inline pga4
pga4_top(void)
{
    return exterior5_top();
}

static inline pga4
pga4_join(pga4 a, pga4 b)
{
    return exterior5_join(a, b);
}

static inline pga4
pga4_meet(pga4 a, pga4 b)
{
    return exterior5_meet(a, b);
}

static inline pga4
pga4_ldual(pga4 a)
{
    return exterior5_ldual(a);
}

static inline pga4
pga4_rdual(pga4 a)
{
    return exterior5_rdual(a);
}

static inline pga4
pga4_dual(pga4 a)
{
    return exterior5_dual(a);
}

/* Central part x + x01234*w of a, as a dual number (w*w = 0). */
static inline dual1
pga4_center_part(pga4 a)
{
    return dual1_make((const double[]){a.data[0], a.data[31]});
}

/*
 * Minimal polynomial over the centre span{1, w}, w = e01234 with w*w = 0,
 * which is the dual1 algebra. Coefficients are in ascending order, p[i]
 * standing for p[i][0] + p[i][1]*w. Write a = c + b with c central and b
 * without central part; b satisfies the centered quartic
 * b^4 - 2*s2*b^2 - (4/3)*s3*b + 2*s2^2 - s4 with central traces s2, s3, s4 of
 * b^2, b^3, b^4, and substituting b = a - c gives the quartic in a.
 */
static inline void
pga4_center_minpoly(pga4 a, dual1 p[5])
{
    dual1 c = pga4_center_part(a), c2 = dual1_square(c), q[3];
    pga4 b = a;
    b.data[0] = b.data[31] = 0;
    pga4 b2 = pga4_square(b), b3 = pga4_mul(b2, b);
    dual1 s2 = pga4_center_part(b2), s3 = pga4_center_part(b3);
    dual1 s4 = pga4_center_part(pga4_square(b2));
    q[0] = vec2_sub(vec2_scale(dual1_square(s2), 2), s4);
    q[1] = vec2_scale(s3, -4./3);
    q[2] = vec2_scale(s2, -2);
    p[0] = vec2_add(vec2_sub(dual1_square(c2), dual1_mul(c, q[1])),
        vec2_add(dual1_mul(c2, q[2]), q[0]));
    p[1] = vec2_sub(vec2_sub(q[1], dual1_mul(vec2_scale(c, 2), q[2])),
        dual1_mul(vec2_scale(c, 4), c2));
    p[2] = vec2_add(vec2_scale(c2, 6), q[2]);
    p[3] = vec2_scale(c, -4);
    p[4] = dual1_one();
}

#endif /* ALGPY_PGA4_H */
