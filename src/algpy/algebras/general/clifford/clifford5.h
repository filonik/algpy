#ifndef ALGPY_CLIFFORD5_H
#define ALGPY_CLIFFORD5_H

#include <string.h>
#include "clifford4.h"
#include "../../../polynomials/hermite_complex.h"
#include "../../../polynomials/hermite8.h"
#include "../../../vectorspaces/vec5.h"
#include "../../../vectorspaces/vec32.h"
#include "clifford1.h"

#define ALGPY_CLIFFORD5_N 32

typedef vec32 clifford5;

static inline clifford5
clifford5_make(const double coeffs[ALGPY_CLIFFORD5_N])
{
    clifford5 result;
    memcpy(result.data, coeffs, sizeof(result.data));
    return result;
}

static inline clifford5
clifford5_grade0(clifford5 a)
{
    const double *x = a.data;
    return clifford5_make((const double[]){x[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline clifford5
clifford5_grade1(clifford5 a)
{
    const double *x = a.data;
    return clifford5_make((const double[]){0, x[1], x[2], 0, x[4], 0, 0, 0, x[8], 0, 0, 0, 0, 0, 0, 0, x[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
}

static inline clifford5
clifford5_grade2(clifford5 a)
{
    const double *x = a.data;
    return clifford5_make((const double[]){0, 0, 0, x[3], 0, x[5], x[6], 0, 0, x[9], x[10], 0, x[12], 0, 0, 0, 0, x[17], x[18], 0, x[20], 0, 0, 0, x[24], 0, 0, 0, 0, 0, 0, 0});
}

static inline clifford5
clifford5_grade3(clifford5 a)
{
    const double *x = a.data;
    return clifford5_make((const double[]){0, 0, 0, 0, 0, 0, 0, x[7], 0, 0, 0, x[11], 0, x[13], x[14], 0, 0, 0, 0, x[19], 0, x[21], x[22], 0, 0, x[25], x[26], 0, x[28], 0, 0, 0});
}

static inline clifford5
clifford5_grade4(clifford5 a)
{
    const double *x = a.data;
    return clifford5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[15], 0, 0, 0, 0, 0, 0, 0, x[23], 0, 0, 0, x[27], 0, x[29], x[30], 0});
}

static inline clifford5
clifford5_grade5(clifford5 a)
{
    const double *x = a.data;
    return clifford5_make((const double[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[31]});
}

static inline clifford5
clifford5_from_re(double a)
{
    clifford5 result = {{0}};
    result.data[0] = a;
    return result;
}

static inline double
clifford5_to_re(clifford5 a)
{
    return a.data[0];
}

static inline clifford5
clifford5_zero(void)
{
    return vec32_zero();
}

static inline clifford5
clifford5_one(void)
{
    return clifford5_from_re(1);
}

static inline clifford5
clifford5_pos(clifford5 a)
{
    return vec32_pos(a);
}

static inline clifford5
clifford5_neg(clifford5 a)
{
    return vec32_neg(a);
}

static inline clifford5
clifford5_add(clifford5 a, clifford5 b)
{
    return vec32_add(a, b);
}

static inline clifford5
clifford5_sub(clifford5 a, clifford5 b)
{
    return vec32_sub(a, b);
}

static inline int
clifford5_equal(clifford5 a, clifford5 b)
{
    return vec32_all_eq(a, b);
}

static inline int
clifford5_not_equal(clifford5 a, clifford5 b)
{
    return !clifford5_equal(a, b);
}

static inline int
clifford5_is_real(clifford5 a)
{
    return vec32_is_real(a);
}

static inline clifford5
clifford5_re(clifford5 a)
{
    return clifford5_from_re(a.data[0]);
}

static inline clifford5
clifford5_im(clifford5 a)
{
    return clifford5_sub(a, clifford5_re(a));
}

static inline clifford5
clifford5_conjugate(clifford5 a)
{
    const double *x = a.data;
    return clifford5_make((const double[]){x[0], -x[1], -x[2], -x[3], -x[4], -x[5], -x[6], x[7], -x[8], -x[9], -x[10], x[11], -x[12], x[13], x[14], x[15], -x[16], -x[17], -x[18], x[19], -x[20], x[21], x[22], x[23], -x[24], x[25], x[26], x[27], x[28], x[29], x[30], -x[31]});
}

static inline void
clifford5_basis(clifford5 result[32])
{
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) result[i].data[j] = i == j ? 1 : 0;
    }
}

static inline void
clifford5_gens(clifford5 result[5])
{
    result[0] = clifford5_make((const double[ALGPY_CLIFFORD5_N]){[1] = 1});
    result[1] = clifford5_make((const double[ALGPY_CLIFFORD5_N]){[2] = 1});
    result[2] = clifford5_make((const double[ALGPY_CLIFFORD5_N]){[4] = 1});
    result[3] = clifford5_make((const double[ALGPY_CLIFFORD5_N]){[8] = 1});
    result[4] = clifford5_make((const double[ALGPY_CLIFFORD5_N]){[16] = 1});
}

static inline clifford5
clifford5_nan(void)
{
    clifford5 result;
    for (int i = 0; i < 32; i++) result.data[i] = NAN;
    return result;
}

static inline clifford5
clifford5_mul(vec5 k, clifford5 a, clifford5 b)
{
    const double *x = a.data, *y = b.data;
    double k0 = k.data[0];
    double k1 = k.data[1];
    double k01 = k0*k1;
    double k2 = k.data[2];
    double k02 = k0*k2;
    double k12 = k1*k2;
    double k012 = k01*k2;
    double k3 = k.data[3];
    double k03 = k0*k3;
    double k13 = k1*k3;
    double k013 = k01*k3;
    double k23 = k2*k3;
    double k023 = k02*k3;
    double k123 = k12*k3;
    double k0123 = k012*k3;
    double k4 = k.data[4];
    double k04 = k0*k4;
    double k14 = k1*k4;
    double k014 = k01*k4;
    double k24 = k2*k4;
    double k024 = k02*k4;
    double k124 = k12*k4;
    double k0124 = k012*k4;
    double k34 = k3*k4;
    double k034 = k03*k4;
    double k134 = k13*k4;
    double k0134 = k013*k4;
    double k234 = k23*k4;
    double k0234 = k023*k4;
    double k1234 = k123*k4;
    double k01234 = k0123*k4;
    return clifford5_make((const double[]){
        x[0]*y[0] + k0*x[1]*y[1] + k1*x[2]*y[2] - k01*x[3]*y[3] + k2*x[4]*y[4] - k02*x[5]*y[5] - k12*x[6]*y[6] - k012*x[7]*y[7] + k3*x[8]*y[8] - k03*x[9]*y[9] - k13*x[10]*y[10] - k013*x[11]*y[11] - k23*x[12]*y[12] - k023*x[13]*y[13] - k123*x[14]*y[14] + k0123*x[15]*y[15] + k4*x[16]*y[16] - k04*x[17]*y[17] - k14*x[18]*y[18] - k014*x[19]*y[19] - k24*x[20]*y[20] - k024*x[21]*y[21] - k124*x[22]*y[22] + k0124*x[23]*y[23] - k34*x[24]*y[24] - k034*x[25]*y[25] - k134*x[26]*y[26] + k0134*x[27]*y[27] - k234*x[28]*y[28] + k0234*x[29]*y[29] + k1234*x[30]*y[30] + k01234*x[31]*y[31],
        x[0]*y[1] + x[1]*y[0] - k1*x[2]*y[3] + k1*x[3]*y[2] - k2*x[4]*y[5] + k2*x[5]*y[4] - k12*x[6]*y[7] - k12*x[7]*y[6] - k3*x[8]*y[9] + k3*x[9]*y[8] - k13*x[10]*y[11] - k13*x[11]*y[10] - k23*x[12]*y[13] - k23*x[13]*y[12] + k123*x[14]*y[15] - k123*x[15]*y[14] - k4*x[16]*y[17] + k4*x[17]*y[16] - k14*x[18]*y[19] - k14*x[19]*y[18] - k24*x[20]*y[21] - k24*x[21]*y[20] + k124*x[22]*y[23] - k124*x[23]*y[22] - k34*x[24]*y[25] - k34*x[25]*y[24] + k134*x[26]*y[27] - k134*x[27]*y[26] + k234*x[28]*y[29] - k234*x[29]*y[28] + k1234*x[30]*y[31] + k1234*x[31]*y[30],
        x[0]*y[2] + k0*x[1]*y[3] + x[2]*y[0] - k0*x[3]*y[1] - k2*x[4]*y[6] + k02*x[5]*y[7] + k2*x[6]*y[4] + k02*x[7]*y[5] - k3*x[8]*y[10] + k03*x[9]*y[11] + k3*x[10]*y[8] + k03*x[11]*y[9] - k23*x[12]*y[14] - k023*x[13]*y[15] - k23*x[14]*y[12] + k023*x[15]*y[13] - k4*x[16]*y[18] + k04*x[17]*y[19] + k4*x[18]*y[16] + k04*x[19]*y[17] - k24*x[20]*y[22] - k024*x[21]*y[23] - k24*x[22]*y[20] + k024*x[23]*y[21] - k34*x[24]*y[26] - k034*x[25]*y[27] - k34*x[26]*y[24] + k034*x[27]*y[25] + k234*x[28]*y[30] - k0234*x[29]*y[31] - k234*x[30]*y[28] - k0234*x[31]*y[29],
        x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0] + k2*x[4]*y[7] - k2*x[5]*y[6] + k2*x[6]*y[5] + k2*x[7]*y[4] + k3*x[8]*y[11] - k3*x[9]*y[10] + k3*x[10]*y[9] + k3*x[11]*y[8] - k23*x[12]*y[15] - k23*x[13]*y[14] + k23*x[14]*y[13] - k23*x[15]*y[12] + k4*x[16]*y[19] - k4*x[17]*y[18] + k4*x[18]*y[17] + k4*x[19]*y[16] - k24*x[20]*y[23] - k24*x[21]*y[22] + k24*x[22]*y[21] - k24*x[23]*y[20] - k34*x[24]*y[27] - k34*x[25]*y[26] + k34*x[26]*y[25] - k34*x[27]*y[24] - k234*x[28]*y[31] + k234*x[29]*y[30] - k234*x[30]*y[29] - k234*x[31]*y[28],
        x[0]*y[4] + k0*x[1]*y[5] + k1*x[2]*y[6] - k01*x[3]*y[7] + x[4]*y[0] - k0*x[5]*y[1] - k1*x[6]*y[2] - k01*x[7]*y[3] - k3*x[8]*y[12] + k03*x[9]*y[13] + k13*x[10]*y[14] + k013*x[11]*y[15] + k3*x[12]*y[8] + k03*x[13]*y[9] + k13*x[14]*y[10] - k013*x[15]*y[11] - k4*x[16]*y[20] + k04*x[17]*y[21] + k14*x[18]*y[22] + k014*x[19]*y[23] + k4*x[20]*y[16] + k04*x[21]*y[17] + k14*x[22]*y[18] - k014*x[23]*y[19] - k34*x[24]*y[28] - k034*x[25]*y[29] - k134*x[26]*y[30] + k0134*x[27]*y[31] - k34*x[28]*y[24] + k034*x[29]*y[25] + k134*x[30]*y[26] + k0134*x[31]*y[27],
        x[0]*y[5] + x[1]*y[4] - k1*x[2]*y[7] + k1*x[3]*y[6] - x[4]*y[1] + x[5]*y[0] - k1*x[6]*y[3] - k1*x[7]*y[2] + k3*x[8]*y[13] - k3*x[9]*y[12] + k13*x[10]*y[15] + k13*x[11]*y[14] + k3*x[12]*y[9] + k3*x[13]*y[8] - k13*x[14]*y[11] + k13*x[15]*y[10] + k4*x[16]*y[21] - k4*x[17]*y[20] + k14*x[18]*y[23] + k14*x[19]*y[22] + k4*x[20]*y[17] + k4*x[21]*y[16] - k14*x[22]*y[19] + k14*x[23]*y[18] - k34*x[24]*y[29] - k34*x[25]*y[28] + k134*x[26]*y[31] - k134*x[27]*y[30] + k34*x[28]*y[25] - k34*x[29]*y[24] + k134*x[30]*y[27] + k134*x[31]*y[26],
        x[0]*y[6] + k0*x[1]*y[7] + x[2]*y[4] - k0*x[3]*y[5] - x[4]*y[2] + k0*x[5]*y[3] + x[6]*y[0] + k0*x[7]*y[1] + k3*x[8]*y[14] - k03*x[9]*y[15] - k3*x[10]*y[12] - k03*x[11]*y[13] + k3*x[12]*y[10] + k03*x[13]*y[11] + k3*x[14]*y[8] - k03*x[15]*y[9] + k4*x[16]*y[22] - k04*x[17]*y[23] - k4*x[18]*y[20] - k04*x[19]*y[21] + k4*x[20]*y[18] + k04*x[21]*y[19] + k4*x[22]*y[16] - k04*x[23]*y[17] - k34*x[24]*y[30] - k034*x[25]*y[31] - k34*x[26]*y[28] + k034*x[27]*y[29] + k34*x[28]*y[26] - k034*x[29]*y[27] - k34*x[30]*y[24] - k034*x[31]*y[25],
        x[0]*y[7] + x[1]*y[6] - x[2]*y[5] + x[3]*y[4] + x[4]*y[3] - x[5]*y[2] + x[6]*y[1] + x[7]*y[0] - k3*x[8]*y[15] + k3*x[9]*y[14] - k3*x[10]*y[13] - k3*x[11]*y[12] + k3*x[12]*y[11] + k3*x[13]*y[10] - k3*x[14]*y[9] + k3*x[15]*y[8] - k4*x[16]*y[23] + k4*x[17]*y[22] - k4*x[18]*y[21] - k4*x[19]*y[20] + k4*x[20]*y[19] + k4*x[21]*y[18] - k4*x[22]*y[17] + k4*x[23]*y[16] - k34*x[24]*y[31] - k34*x[25]*y[30] + k34*x[26]*y[29] - k34*x[27]*y[28] - k34*x[28]*y[27] + k34*x[29]*y[26] - k34*x[30]*y[25] - k34*x[31]*y[24],
        x[0]*y[8] + k0*x[1]*y[9] + k1*x[2]*y[10] - k01*x[3]*y[11] + k2*x[4]*y[12] - k02*x[5]*y[13] - k12*x[6]*y[14] - k012*x[7]*y[15] + x[8]*y[0] - k0*x[9]*y[1] - k1*x[10]*y[2] - k01*x[11]*y[3] - k2*x[12]*y[4] - k02*x[13]*y[5] - k12*x[14]*y[6] + k012*x[15]*y[7] - k4*x[16]*y[24] + k04*x[17]*y[25] + k14*x[18]*y[26] + k014*x[19]*y[27] + k24*x[20]*y[28] + k024*x[21]*y[29] + k124*x[22]*y[30] - k0124*x[23]*y[31] + k4*x[24]*y[16] + k04*x[25]*y[17] + k14*x[26]*y[18] - k014*x[27]*y[19] + k24*x[28]*y[20] - k024*x[29]*y[21] - k124*x[30]*y[22] - k0124*x[31]*y[23],
        x[0]*y[9] + x[1]*y[8] - k1*x[2]*y[11] + k1*x[3]*y[10] - k2*x[4]*y[13] + k2*x[5]*y[12] - k12*x[6]*y[15] - k12*x[7]*y[14] - x[8]*y[1] + x[9]*y[0] - k1*x[10]*y[3] - k1*x[11]*y[2] - k2*x[12]*y[5] - k2*x[13]*y[4] + k12*x[14]*y[7] - k12*x[15]*y[6] + k4*x[16]*y[25] - k4*x[17]*y[24] + k14*x[18]*y[27] + k14*x[19]*y[26] + k24*x[20]*y[29] + k24*x[21]*y[28] - k124*x[22]*y[31] + k124*x[23]*y[30] + k4*x[24]*y[17] + k4*x[25]*y[16] - k14*x[26]*y[19] + k14*x[27]*y[18] - k24*x[28]*y[21] + k24*x[29]*y[20] - k124*x[30]*y[23] - k124*x[31]*y[22],
        x[0]*y[10] + k0*x[1]*y[11] + x[2]*y[8] - k0*x[3]*y[9] - k2*x[4]*y[14] + k02*x[5]*y[15] + k2*x[6]*y[12] + k02*x[7]*y[13] - x[8]*y[2] + k0*x[9]*y[3] + x[10]*y[0] + k0*x[11]*y[1] - k2*x[12]*y[6] - k02*x[13]*y[7] - k2*x[14]*y[4] + k02*x[15]*y[5] + k4*x[16]*y[26] - k04*x[17]*y[27] - k4*x[18]*y[24] - k04*x[19]*y[25] + k24*x[20]*y[30] + k024*x[21]*y[31] + k24*x[22]*y[28] - k024*x[23]*y[29] + k4*x[24]*y[18] + k04*x[25]*y[19] + k4*x[26]*y[16] - k04*x[27]*y[17] - k24*x[28]*y[22] + k024*x[29]*y[23] + k24*x[30]*y[20] + k024*x[31]*y[21],
        x[0]*y[11] + x[1]*y[10] - x[2]*y[9] + x[3]*y[8] + k2*x[4]*y[15] - k2*x[5]*y[14] + k2*x[6]*y[13] + k2*x[7]*y[12] + x[8]*y[3] - x[9]*y[2] + x[10]*y[1] + x[11]*y[0] - k2*x[12]*y[7] - k2*x[13]*y[6] + k2*x[14]*y[5] - k2*x[15]*y[4] - k4*x[16]*y[27] + k4*x[17]*y[26] - k4*x[18]*y[25] - k4*x[19]*y[24] + k24*x[20]*y[31] + k24*x[21]*y[30] - k24*x[22]*y[29] + k24*x[23]*y[28] + k4*x[24]*y[19] + k4*x[25]*y[18] - k4*x[26]*y[17] + k4*x[27]*y[16] + k24*x[28]*y[23] - k24*x[29]*y[22] + k24*x[30]*y[21] + k24*x[31]*y[20],
        x[0]*y[12] + k0*x[1]*y[13] + k1*x[2]*y[14] - k01*x[3]*y[15] + x[4]*y[8] - k0*x[5]*y[9] - k1*x[6]*y[10] - k01*x[7]*y[11] - x[8]*y[4] + k0*x[9]*y[5] + k1*x[10]*y[6] + k01*x[11]*y[7] + x[12]*y[0] + k0*x[13]*y[1] + k1*x[14]*y[2] - k01*x[15]*y[3] + k4*x[16]*y[28] - k04*x[17]*y[29] - k14*x[18]*y[30] - k014*x[19]*y[31] - k4*x[20]*y[24] - k04*x[21]*y[25] - k14*x[22]*y[26] + k014*x[23]*y[27] + k4*x[24]*y[20] + k04*x[25]*y[21] + k14*x[26]*y[22] - k014*x[27]*y[23] + k4*x[28]*y[16] - k04*x[29]*y[17] - k14*x[30]*y[18] - k014*x[31]*y[19],
        x[0]*y[13] + x[1]*y[12] - k1*x[2]*y[15] + k1*x[3]*y[14] - x[4]*y[9] + x[5]*y[8] - k1*x[6]*y[11] - k1*x[7]*y[10] + x[8]*y[5] - x[9]*y[4] + k1*x[10]*y[7] + k1*x[11]*y[6] + x[12]*y[1] + x[13]*y[0] - k1*x[14]*y[3] + k1*x[15]*y[2] - k4*x[16]*y[29] + k4*x[17]*y[28] - k14*x[18]*y[31] - k14*x[19]*y[30] - k4*x[20]*y[25] - k4*x[21]*y[24] + k14*x[22]*y[27] - k14*x[23]*y[26] + k4*x[24]*y[21] + k4*x[25]*y[20] - k14*x[26]*y[23] + k14*x[27]*y[22] - k4*x[28]*y[17] + k4*x[29]*y[16] - k14*x[30]*y[19] - k14*x[31]*y[18],
        x[0]*y[14] + k0*x[1]*y[15] + x[2]*y[12] - k0*x[3]*y[13] - x[4]*y[10] + k0*x[5]*y[11] + x[6]*y[8] + k0*x[7]*y[9] + x[8]*y[6] - k0*x[9]*y[7] - x[10]*y[4] - k0*x[11]*y[5] + x[12]*y[2] + k0*x[13]*y[3] + x[14]*y[0] - k0*x[15]*y[1] - k4*x[16]*y[30] + k04*x[17]*y[31] + k4*x[18]*y[28] + k04*x[19]*y[29] - k4*x[20]*y[26] - k04*x[21]*y[27] - k4*x[22]*y[24] + k04*x[23]*y[25] + k4*x[24]*y[22] + k04*x[25]*y[23] + k4*x[26]*y[20] - k04*x[27]*y[21] - k4*x[28]*y[18] + k04*x[29]*y[19] + k4*x[30]*y[16] + k04*x[31]*y[17],
        x[0]*y[15] + x[1]*y[14] - x[2]*y[13] + x[3]*y[12] + x[4]*y[11] - x[5]*y[10] + x[6]*y[9] + x[7]*y[8] - x[8]*y[7] + x[9]*y[6] - x[10]*y[5] - x[11]*y[4] + x[12]*y[3] + x[13]*y[2] - x[14]*y[1] + x[15]*y[0] + k4*x[16]*y[31] - k4*x[17]*y[30] + k4*x[18]*y[29] + k4*x[19]*y[28] - k4*x[20]*y[27] - k4*x[21]*y[26] + k4*x[22]*y[25] - k4*x[23]*y[24] + k4*x[24]*y[23] + k4*x[25]*y[22] - k4*x[26]*y[21] + k4*x[27]*y[20] + k4*x[28]*y[19] - k4*x[29]*y[18] + k4*x[30]*y[17] + k4*x[31]*y[16],
        x[0]*y[16] + k0*x[1]*y[17] + k1*x[2]*y[18] - k01*x[3]*y[19] + k2*x[4]*y[20] - k02*x[5]*y[21] - k12*x[6]*y[22] - k012*x[7]*y[23] + k3*x[8]*y[24] - k03*x[9]*y[25] - k13*x[10]*y[26] - k013*x[11]*y[27] - k23*x[12]*y[28] - k023*x[13]*y[29] - k123*x[14]*y[30] + k0123*x[15]*y[31] + x[16]*y[0] - k0*x[17]*y[1] - k1*x[18]*y[2] - k01*x[19]*y[3] - k2*x[20]*y[4] - k02*x[21]*y[5] - k12*x[22]*y[6] + k012*x[23]*y[7] - k3*x[24]*y[8] - k03*x[25]*y[9] - k13*x[26]*y[10] + k013*x[27]*y[11] - k23*x[28]*y[12] + k023*x[29]*y[13] + k123*x[30]*y[14] + k0123*x[31]*y[15],
        x[0]*y[17] + x[1]*y[16] - k1*x[2]*y[19] + k1*x[3]*y[18] - k2*x[4]*y[21] + k2*x[5]*y[20] - k12*x[6]*y[23] - k12*x[7]*y[22] - k3*x[8]*y[25] + k3*x[9]*y[24] - k13*x[10]*y[27] - k13*x[11]*y[26] - k23*x[12]*y[29] - k23*x[13]*y[28] + k123*x[14]*y[31] - k123*x[15]*y[30] - x[16]*y[1] + x[17]*y[0] - k1*x[18]*y[3] - k1*x[19]*y[2] - k2*x[20]*y[5] - k2*x[21]*y[4] + k12*x[22]*y[7] - k12*x[23]*y[6] - k3*x[24]*y[9] - k3*x[25]*y[8] + k13*x[26]*y[11] - k13*x[27]*y[10] + k23*x[28]*y[13] - k23*x[29]*y[12] + k123*x[30]*y[15] + k123*x[31]*y[14],
        x[0]*y[18] + k0*x[1]*y[19] + x[2]*y[16] - k0*x[3]*y[17] - k2*x[4]*y[22] + k02*x[5]*y[23] + k2*x[6]*y[20] + k02*x[7]*y[21] - k3*x[8]*y[26] + k03*x[9]*y[27] + k3*x[10]*y[24] + k03*x[11]*y[25] - k23*x[12]*y[30] - k023*x[13]*y[31] - k23*x[14]*y[28] + k023*x[15]*y[29] - x[16]*y[2] + k0*x[17]*y[3] + x[18]*y[0] + k0*x[19]*y[1] - k2*x[20]*y[6] - k02*x[21]*y[7] - k2*x[22]*y[4] + k02*x[23]*y[5] - k3*x[24]*y[10] - k03*x[25]*y[11] - k3*x[26]*y[8] + k03*x[27]*y[9] + k23*x[28]*y[14] - k023*x[29]*y[15] - k23*x[30]*y[12] - k023*x[31]*y[13],
        x[0]*y[19] + x[1]*y[18] - x[2]*y[17] + x[3]*y[16] + k2*x[4]*y[23] - k2*x[5]*y[22] + k2*x[6]*y[21] + k2*x[7]*y[20] + k3*x[8]*y[27] - k3*x[9]*y[26] + k3*x[10]*y[25] + k3*x[11]*y[24] - k23*x[12]*y[31] - k23*x[13]*y[30] + k23*x[14]*y[29] - k23*x[15]*y[28] + x[16]*y[3] - x[17]*y[2] + x[18]*y[1] + x[19]*y[0] - k2*x[20]*y[7] - k2*x[21]*y[6] + k2*x[22]*y[5] - k2*x[23]*y[4] - k3*x[24]*y[11] - k3*x[25]*y[10] + k3*x[26]*y[9] - k3*x[27]*y[8] - k23*x[28]*y[15] + k23*x[29]*y[14] - k23*x[30]*y[13] - k23*x[31]*y[12],
        x[0]*y[20] + k0*x[1]*y[21] + k1*x[2]*y[22] - k01*x[3]*y[23] + x[4]*y[16] - k0*x[5]*y[17] - k1*x[6]*y[18] - k01*x[7]*y[19] - k3*x[8]*y[28] + k03*x[9]*y[29] + k13*x[10]*y[30] + k013*x[11]*y[31] + k3*x[12]*y[24] + k03*x[13]*y[25] + k13*x[14]*y[26] - k013*x[15]*y[27] - x[16]*y[4] + k0*x[17]*y[5] + k1*x[18]*y[6] + k01*x[19]*y[7] + x[20]*y[0] + k0*x[21]*y[1] + k1*x[22]*y[2] - k01*x[23]*y[3] - k3*x[24]*y[12] - k03*x[25]*y[13] - k13*x[26]*y[14] + k013*x[27]*y[15] - k3*x[28]*y[8] + k03*x[29]*y[9] + k13*x[30]*y[10] + k013*x[31]*y[11],
        x[0]*y[21] + x[1]*y[20] - k1*x[2]*y[23] + k1*x[3]*y[22] - x[4]*y[17] + x[5]*y[16] - k1*x[6]*y[19] - k1*x[7]*y[18] + k3*x[8]*y[29] - k3*x[9]*y[28] + k13*x[10]*y[31] + k13*x[11]*y[30] + k3*x[12]*y[25] + k3*x[13]*y[24] - k13*x[14]*y[27] + k13*x[15]*y[26] + x[16]*y[5] - x[17]*y[4] + k1*x[18]*y[7] + k1*x[19]*y[6] + x[20]*y[1] + x[21]*y[0] - k1*x[22]*y[3] + k1*x[23]*y[2] - k3*x[24]*y[13] - k3*x[25]*y[12] + k13*x[26]*y[15] - k13*x[27]*y[14] + k3*x[28]*y[9] - k3*x[29]*y[8] + k13*x[30]*y[11] + k13*x[31]*y[10],
        x[0]*y[22] + k0*x[1]*y[23] + x[2]*y[20] - k0*x[3]*y[21] - x[4]*y[18] + k0*x[5]*y[19] + x[6]*y[16] + k0*x[7]*y[17] + k3*x[8]*y[30] - k03*x[9]*y[31] - k3*x[10]*y[28] - k03*x[11]*y[29] + k3*x[12]*y[26] + k03*x[13]*y[27] + k3*x[14]*y[24] - k03*x[15]*y[25] + x[16]*y[6] - k0*x[17]*y[7] - x[18]*y[4] - k0*x[19]*y[5] + x[20]*y[2] + k0*x[21]*y[3] + x[22]*y[0] - k0*x[23]*y[1] - k3*x[24]*y[14] - k03*x[25]*y[15] - k3*x[26]*y[12] + k03*x[27]*y[13] + k3*x[28]*y[10] - k03*x[29]*y[11] - k3*x[30]*y[8] - k03*x[31]*y[9],
        x[0]*y[23] + x[1]*y[22] - x[2]*y[21] + x[3]*y[20] + x[4]*y[19] - x[5]*y[18] + x[6]*y[17] + x[7]*y[16] - k3*x[8]*y[31] + k3*x[9]*y[30] - k3*x[10]*y[29] - k3*x[11]*y[28] + k3*x[12]*y[27] + k3*x[13]*y[26] - k3*x[14]*y[25] + k3*x[15]*y[24] - x[16]*y[7] + x[17]*y[6] - x[18]*y[5] - x[19]*y[4] + x[20]*y[3] + x[21]*y[2] - x[22]*y[1] + x[23]*y[0] - k3*x[24]*y[15] - k3*x[25]*y[14] + k3*x[26]*y[13] - k3*x[27]*y[12] - k3*x[28]*y[11] + k3*x[29]*y[10] - k3*x[30]*y[9] - k3*x[31]*y[8],
        x[0]*y[24] + k0*x[1]*y[25] + k1*x[2]*y[26] - k01*x[3]*y[27] + k2*x[4]*y[28] - k02*x[5]*y[29] - k12*x[6]*y[30] - k012*x[7]*y[31] + x[8]*y[16] - k0*x[9]*y[17] - k1*x[10]*y[18] - k01*x[11]*y[19] - k2*x[12]*y[20] - k02*x[13]*y[21] - k12*x[14]*y[22] + k012*x[15]*y[23] - x[16]*y[8] + k0*x[17]*y[9] + k1*x[18]*y[10] + k01*x[19]*y[11] + k2*x[20]*y[12] + k02*x[21]*y[13] + k12*x[22]*y[14] - k012*x[23]*y[15] + x[24]*y[0] + k0*x[25]*y[1] + k1*x[26]*y[2] - k01*x[27]*y[3] + k2*x[28]*y[4] - k02*x[29]*y[5] - k12*x[30]*y[6] - k012*x[31]*y[7],
        x[0]*y[25] + x[1]*y[24] - k1*x[2]*y[27] + k1*x[3]*y[26] - k2*x[4]*y[29] + k2*x[5]*y[28] - k12*x[6]*y[31] - k12*x[7]*y[30] - x[8]*y[17] + x[9]*y[16] - k1*x[10]*y[19] - k1*x[11]*y[18] - k2*x[12]*y[21] - k2*x[13]*y[20] + k12*x[14]*y[23] - k12*x[15]*y[22] + x[16]*y[9] - x[17]*y[8] + k1*x[18]*y[11] + k1*x[19]*y[10] + k2*x[20]*y[13] + k2*x[21]*y[12] - k12*x[22]*y[15] + k12*x[23]*y[14] + x[24]*y[1] + x[25]*y[0] - k1*x[26]*y[3] + k1*x[27]*y[2] - k2*x[28]*y[5] + k2*x[29]*y[4] - k12*x[30]*y[7] - k12*x[31]*y[6],
        x[0]*y[26] + k0*x[1]*y[27] + x[2]*y[24] - k0*x[3]*y[25] - k2*x[4]*y[30] + k02*x[5]*y[31] + k2*x[6]*y[28] + k02*x[7]*y[29] - x[8]*y[18] + k0*x[9]*y[19] + x[10]*y[16] + k0*x[11]*y[17] - k2*x[12]*y[22] - k02*x[13]*y[23] - k2*x[14]*y[20] + k02*x[15]*y[21] + x[16]*y[10] - k0*x[17]*y[11] - x[18]*y[8] - k0*x[19]*y[9] + k2*x[20]*y[14] + k02*x[21]*y[15] + k2*x[22]*y[12] - k02*x[23]*y[13] + x[24]*y[2] + k0*x[25]*y[3] + x[26]*y[0] - k0*x[27]*y[1] - k2*x[28]*y[6] + k02*x[29]*y[7] + k2*x[30]*y[4] + k02*x[31]*y[5],
        x[0]*y[27] + x[1]*y[26] - x[2]*y[25] + x[3]*y[24] + k2*x[4]*y[31] - k2*x[5]*y[30] + k2*x[6]*y[29] + k2*x[7]*y[28] + x[8]*y[19] - x[9]*y[18] + x[10]*y[17] + x[11]*y[16] - k2*x[12]*y[23] - k2*x[13]*y[22] + k2*x[14]*y[21] - k2*x[15]*y[20] - x[16]*y[11] + x[17]*y[10] - x[18]*y[9] - x[19]*y[8] + k2*x[20]*y[15] + k2*x[21]*y[14] - k2*x[22]*y[13] + k2*x[23]*y[12] + x[24]*y[3] + x[25]*y[2] - x[26]*y[1] + x[27]*y[0] + k2*x[28]*y[7] - k2*x[29]*y[6] + k2*x[30]*y[5] + k2*x[31]*y[4],
        x[0]*y[28] + k0*x[1]*y[29] + k1*x[2]*y[30] - k01*x[3]*y[31] + x[4]*y[24] - k0*x[5]*y[25] - k1*x[6]*y[26] - k01*x[7]*y[27] - x[8]*y[20] + k0*x[9]*y[21] + k1*x[10]*y[22] + k01*x[11]*y[23] + x[12]*y[16] + k0*x[13]*y[17] + k1*x[14]*y[18] - k01*x[15]*y[19] + x[16]*y[12] - k0*x[17]*y[13] - k1*x[18]*y[14] - k01*x[19]*y[15] - x[20]*y[8] - k0*x[21]*y[9] - k1*x[22]*y[10] + k01*x[23]*y[11] + x[24]*y[4] + k0*x[25]*y[5] + k1*x[26]*y[6] - k01*x[27]*y[7] + x[28]*y[0] - k0*x[29]*y[1] - k1*x[30]*y[2] - k01*x[31]*y[3],
        x[0]*y[29] + x[1]*y[28] - k1*x[2]*y[31] + k1*x[3]*y[30] - x[4]*y[25] + x[5]*y[24] - k1*x[6]*y[27] - k1*x[7]*y[26] + x[8]*y[21] - x[9]*y[20] + k1*x[10]*y[23] + k1*x[11]*y[22] + x[12]*y[17] + x[13]*y[16] - k1*x[14]*y[19] + k1*x[15]*y[18] - x[16]*y[13] + x[17]*y[12] - k1*x[18]*y[15] - k1*x[19]*y[14] - x[20]*y[9] - x[21]*y[8] + k1*x[22]*y[11] - k1*x[23]*y[10] + x[24]*y[5] + x[25]*y[4] - k1*x[26]*y[7] + k1*x[27]*y[6] - x[28]*y[1] + x[29]*y[0] - k1*x[30]*y[3] - k1*x[31]*y[2],
        x[0]*y[30] + k0*x[1]*y[31] + x[2]*y[28] - k0*x[3]*y[29] - x[4]*y[26] + k0*x[5]*y[27] + x[6]*y[24] + k0*x[7]*y[25] + x[8]*y[22] - k0*x[9]*y[23] - x[10]*y[20] - k0*x[11]*y[21] + x[12]*y[18] + k0*x[13]*y[19] + x[14]*y[16] - k0*x[15]*y[17] - x[16]*y[14] + k0*x[17]*y[15] + x[18]*y[12] + k0*x[19]*y[13] - x[20]*y[10] - k0*x[21]*y[11] - x[22]*y[8] + k0*x[23]*y[9] + x[24]*y[6] + k0*x[25]*y[7] + x[26]*y[4] - k0*x[27]*y[5] - x[28]*y[2] + k0*x[29]*y[3] + x[30]*y[0] + k0*x[31]*y[1],
        x[0]*y[31] + x[1]*y[30] - x[2]*y[29] + x[3]*y[28] + x[4]*y[27] - x[5]*y[26] + x[6]*y[25] + x[7]*y[24] - x[8]*y[23] + x[9]*y[22] - x[10]*y[21] - x[11]*y[20] + x[12]*y[19] + x[13]*y[18] - x[14]*y[17] + x[15]*y[16] + x[16]*y[15] - x[17]*y[14] + x[18]*y[13] + x[19]*y[12] - x[20]*y[11] - x[21]*y[10] + x[22]*y[9] - x[23]*y[8] + x[24]*y[7] + x[25]*y[6] - x[26]*y[5] + x[27]*y[4] + x[28]*y[3] - x[29]*y[2] + x[30]*y[1] + x[31]*y[0]});
}

static inline clifford5
clifford5_square(vec5 k, clifford5 a)
{
    const double *x = a.data;
    double k0 = k.data[0];
    double k1 = k.data[1];
    double k01 = k0*k1;
    double k2 = k.data[2];
    double k02 = k0*k2;
    double k12 = k1*k2;
    double k012 = k01*k2;
    double k3 = k.data[3];
    double k03 = k0*k3;
    double k13 = k1*k3;
    double k013 = k01*k3;
    double k23 = k2*k3;
    double k023 = k02*k3;
    double k123 = k12*k3;
    double k0123 = k012*k3;
    double k4 = k.data[4];
    double k04 = k0*k4;
    double k14 = k1*k4;
    double k014 = k01*k4;
    double k24 = k2*k4;
    double k024 = k02*k4;
    double k124 = k12*k4;
    double k0124 = k012*k4;
    double k34 = k3*k4;
    double k034 = k03*k4;
    double k134 = k13*k4;
    double k0134 = k013*k4;
    double k234 = k23*k4;
    double k0234 = k023*k4;
    double k1234 = k123*k4;
    double k01234 = k0123*k4;
    return clifford5_make((const double[]){
        x[0]*x[0] + k0*x[1]*x[1] + k1*x[2]*x[2] - k01*x[3]*x[3] + k2*x[4]*x[4] - k02*x[5]*x[5] - k12*x[6]*x[6] - k012*x[7]*x[7] + k3*x[8]*x[8] - k03*x[9]*x[9] - k13*x[10]*x[10] - k013*x[11]*x[11] - k23*x[12]*x[12] - k023*x[13]*x[13] - k123*x[14]*x[14] + k0123*x[15]*x[15] + k4*x[16]*x[16] - k04*x[17]*x[17] - k14*x[18]*x[18] - k014*x[19]*x[19] - k24*x[20]*x[20] - k024*x[21]*x[21] - k124*x[22]*x[22] + k0124*x[23]*x[23] - k34*x[24]*x[24] - k034*x[25]*x[25] - k134*x[26]*x[26] + k0134*x[27]*x[27] - k234*x[28]*x[28] + k0234*x[29]*x[29] + k1234*x[30]*x[30] + k01234*x[31]*x[31],
        2*x[0]*x[1] - 2*k12*x[6]*x[7] - 2*k13*x[10]*x[11] - 2*k23*x[12]*x[13] - 2*k14*x[18]*x[19] - 2*k24*x[20]*x[21] - 2*k34*x[24]*x[25] + 2*k1234*x[30]*x[31],
        2*x[0]*x[2] + 2*k02*x[5]*x[7] + 2*k03*x[9]*x[11] - 2*k23*x[12]*x[14] + 2*k04*x[17]*x[19] - 2*k24*x[20]*x[22] - 2*k34*x[24]*x[26] - 2*k0234*x[29]*x[31],
        2*x[0]*x[3] + 2*k2*x[4]*x[7] + 2*k3*x[8]*x[11] - 2*k23*x[12]*x[15] + 2*k4*x[16]*x[19] - 2*k24*x[20]*x[23] - 2*k34*x[24]*x[27] - 2*k234*x[28]*x[31],
        2*x[0]*x[4] - 2*k01*x[3]*x[7] + 2*k03*x[9]*x[13] + 2*k13*x[10]*x[14] + 2*k04*x[17]*x[21] + 2*k14*x[18]*x[22] - 2*k34*x[24]*x[28] + 2*k0134*x[27]*x[31],
        2*x[0]*x[5] - 2*k1*x[2]*x[7] + 2*k3*x[8]*x[13] + 2*k13*x[10]*x[15] + 2*k4*x[16]*x[21] + 2*k14*x[18]*x[23] - 2*k34*x[24]*x[29] + 2*k134*x[26]*x[31],
        2*x[0]*x[6] + 2*k0*x[1]*x[7] + 2*k3*x[8]*x[14] - 2*k03*x[9]*x[15] + 2*k4*x[16]*x[22] - 2*k04*x[17]*x[23] - 2*k34*x[24]*x[30] - 2*k034*x[25]*x[31],
        2*x[0]*x[7] + 2*x[1]*x[6] - 2*x[2]*x[5] + 2*x[3]*x[4] - 2*k34*x[24]*x[31] - 2*k34*x[25]*x[30] + 2*k34*x[26]*x[29] - 2*k34*x[27]*x[28],
        2*x[0]*x[8] - 2*k01*x[3]*x[11] - 2*k02*x[5]*x[13] - 2*k12*x[6]*x[14] + 2*k04*x[17]*x[25] + 2*k14*x[18]*x[26] + 2*k24*x[20]*x[28] - 2*k0124*x[23]*x[31],
        2*x[0]*x[9] - 2*k1*x[2]*x[11] - 2*k2*x[4]*x[13] - 2*k12*x[6]*x[15] + 2*k4*x[16]*x[25] + 2*k14*x[18]*x[27] + 2*k24*x[20]*x[29] - 2*k124*x[22]*x[31],
        2*x[0]*x[10] + 2*k0*x[1]*x[11] - 2*k2*x[4]*x[14] + 2*k02*x[5]*x[15] + 2*k4*x[16]*x[26] - 2*k04*x[17]*x[27] + 2*k24*x[20]*x[30] + 2*k024*x[21]*x[31],
        2*x[0]*x[11] + 2*x[1]*x[10] - 2*x[2]*x[9] + 2*x[3]*x[8] + 2*k24*x[20]*x[31] + 2*k24*x[21]*x[30] - 2*k24*x[22]*x[29] + 2*k24*x[23]*x[28],
        2*x[0]*x[12] + 2*k0*x[1]*x[13] + 2*k1*x[2]*x[14] - 2*k01*x[3]*x[15] + 2*k4*x[16]*x[28] - 2*k04*x[17]*x[29] - 2*k14*x[18]*x[30] - 2*k014*x[19]*x[31],
        2*x[0]*x[13] + 2*x[1]*x[12] - 2*x[4]*x[9] + 2*x[5]*x[8] - 2*k14*x[18]*x[31] - 2*k14*x[19]*x[30] + 2*k14*x[22]*x[27] - 2*k14*x[23]*x[26],
        2*x[0]*x[14] + 2*x[2]*x[12] - 2*x[4]*x[10] + 2*x[6]*x[8] + 2*k04*x[17]*x[31] + 2*k04*x[19]*x[29] - 2*k04*x[21]*x[27] + 2*k04*x[23]*x[25],
        2*x[0]*x[15] + 2*x[3]*x[12] - 2*x[5]*x[10] + 2*x[6]*x[9] + 2*k4*x[16]*x[31] + 2*k4*x[19]*x[28] - 2*k4*x[21]*x[26] + 2*k4*x[22]*x[25],
        2*x[0]*x[16] - 2*k01*x[3]*x[19] - 2*k02*x[5]*x[21] - 2*k12*x[6]*x[22] - 2*k03*x[9]*x[25] - 2*k13*x[10]*x[26] - 2*k23*x[12]*x[28] + 2*k0123*x[15]*x[31],
        2*x[0]*x[17] - 2*k1*x[2]*x[19] - 2*k2*x[4]*x[21] - 2*k12*x[6]*x[23] - 2*k3*x[8]*x[25] - 2*k13*x[10]*x[27] - 2*k23*x[12]*x[29] + 2*k123*x[14]*x[31],
        2*x[0]*x[18] + 2*k0*x[1]*x[19] - 2*k2*x[4]*x[22] + 2*k02*x[5]*x[23] - 2*k3*x[8]*x[26] + 2*k03*x[9]*x[27] - 2*k23*x[12]*x[30] - 2*k023*x[13]*x[31],
        2*x[0]*x[19] + 2*x[1]*x[18] - 2*x[2]*x[17] + 2*x[3]*x[16] - 2*k23*x[12]*x[31] - 2*k23*x[13]*x[30] + 2*k23*x[14]*x[29] - 2*k23*x[15]*x[28],
        2*x[0]*x[20] + 2*k0*x[1]*x[21] + 2*k1*x[2]*x[22] - 2*k01*x[3]*x[23] - 2*k3*x[8]*x[28] + 2*k03*x[9]*x[29] + 2*k13*x[10]*x[30] + 2*k013*x[11]*x[31],
        2*x[0]*x[21] + 2*x[1]*x[20] - 2*x[4]*x[17] + 2*x[5]*x[16] + 2*k13*x[10]*x[31] + 2*k13*x[11]*x[30] - 2*k13*x[14]*x[27] + 2*k13*x[15]*x[26],
        2*x[0]*x[22] + 2*x[2]*x[20] - 2*x[4]*x[18] + 2*x[6]*x[16] - 2*k03*x[9]*x[31] - 2*k03*x[11]*x[29] + 2*k03*x[13]*x[27] - 2*k03*x[15]*x[25],
        2*x[0]*x[23] + 2*x[3]*x[20] - 2*x[5]*x[18] + 2*x[6]*x[17] - 2*k3*x[8]*x[31] - 2*k3*x[11]*x[28] + 2*k3*x[13]*x[26] - 2*k3*x[14]*x[25],
        2*x[0]*x[24] + 2*k0*x[1]*x[25] + 2*k1*x[2]*x[26] - 2*k01*x[3]*x[27] + 2*k2*x[4]*x[28] - 2*k02*x[5]*x[29] - 2*k12*x[6]*x[30] - 2*k012*x[7]*x[31],
        2*x[0]*x[25] + 2*x[1]*x[24] - 2*k12*x[6]*x[31] - 2*k12*x[7]*x[30] - 2*x[8]*x[17] + 2*x[9]*x[16] + 2*k12*x[14]*x[23] - 2*k12*x[15]*x[22],
        2*x[0]*x[26] + 2*x[2]*x[24] + 2*k02*x[5]*x[31] + 2*k02*x[7]*x[29] - 2*x[8]*x[18] + 2*x[10]*x[16] - 2*k02*x[13]*x[23] + 2*k02*x[15]*x[21],
        2*x[0]*x[27] + 2*x[3]*x[24] + 2*k2*x[4]*x[31] + 2*k2*x[7]*x[28] - 2*x[9]*x[18] + 2*x[10]*x[17] - 2*k2*x[13]*x[22] + 2*k2*x[14]*x[21],
        2*x[0]*x[28] - 2*k01*x[3]*x[31] + 2*x[4]*x[24] - 2*k01*x[7]*x[27] - 2*x[8]*x[20] + 2*k01*x[11]*x[23] + 2*x[12]*x[16] - 2*k01*x[15]*x[19],
        2*x[0]*x[29] - 2*k1*x[2]*x[31] + 2*x[5]*x[24] - 2*k1*x[7]*x[26] - 2*x[9]*x[20] + 2*k1*x[11]*x[22] + 2*x[12]*x[17] - 2*k1*x[14]*x[19],
        2*x[0]*x[30] + 2*k0*x[1]*x[31] + 2*x[6]*x[24] + 2*k0*x[7]*x[25] - 2*x[10]*x[20] - 2*k0*x[11]*x[21] + 2*x[12]*x[18] + 2*k0*x[13]*x[19],
        2*x[0]*x[31] + 2*x[1]*x[30] - 2*x[2]*x[29] + 2*x[3]*x[28] + 2*x[4]*x[27] - 2*x[5]*x[26] + 2*x[6]*x[25] + 2*x[7]*x[24] - 2*x[8]*x[23] + 2*x[9]*x[22] - 2*x[10]*x[21] - 2*x[11]*x[20] + 2*x[12]*x[19] + 2*x[13]*x[18] - 2*x[14]*x[17] + 2*x[15]*x[16]});
}

static inline void
clifford5_coeffs(vec5 k, clifford5 result[1024])
{
    clifford5 basis[32];
    clifford5_basis(basis);
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) result[i*32+j] = clifford5_mul(k, basis[i], basis[j]);
    }
}

static inline double
clifford5_scalar_product(vec5 k, clifford5 a, clifford5 b)
{
    const double *x = a.data, *y = b.data;
    double k0 = k.data[0];
    double k1 = k.data[1];
    double k01 = k0*k1;
    double k2 = k.data[2];
    double k02 = k0*k2;
    double k12 = k1*k2;
    double k012 = k01*k2;
    double k3 = k.data[3];
    double k03 = k0*k3;
    double k13 = k1*k3;
    double k013 = k01*k3;
    double k23 = k2*k3;
    double k023 = k02*k3;
    double k123 = k12*k3;
    double k0123 = k012*k3;
    double k4 = k.data[4];
    double k04 = k0*k4;
    double k14 = k1*k4;
    double k014 = k01*k4;
    double k24 = k2*k4;
    double k024 = k02*k4;
    double k124 = k12*k4;
    double k0124 = k012*k4;
    double k34 = k3*k4;
    double k034 = k03*k4;
    double k134 = k13*k4;
    double k0134 = k013*k4;
    double k234 = k23*k4;
    double k0234 = k023*k4;
    double k1234 = k123*k4;
    double k01234 = k0123*k4;
    return x[0]*y[0] + k0*x[1]*y[1] + k1*x[2]*y[2] - k01*x[3]*y[3] + k2*x[4]*y[4] - k02*x[5]*y[5] - k12*x[6]*y[6] - k012*x[7]*y[7] + k3*x[8]*y[8] - k03*x[9]*y[9] - k13*x[10]*y[10] - k013*x[11]*y[11] - k23*x[12]*y[12] - k023*x[13]*y[13] - k123*x[14]*y[14] + k0123*x[15]*y[15] + k4*x[16]*y[16] - k04*x[17]*y[17] - k14*x[18]*y[18] - k014*x[19]*y[19] - k24*x[20]*y[20] - k024*x[21]*y[21] - k124*x[22]*y[22] + k0124*x[23]*y[23] - k34*x[24]*y[24] - k034*x[25]*y[25] - k134*x[26]*y[26] + k0134*x[27]*y[27] - k234*x[28]*y[28] + k0234*x[29]*y[29] + k1234*x[30]*y[30] + k01234*x[31]*y[31];
}

static inline double
clifford5_pseudoscalar_square(vec5 k)
{
    return k.data[0]*k.data[1]*k.data[2]*k.data[3]*k.data[4];
}

static inline double
clifford5_pseudoscalar_product(vec5 k, int mask)
{
    static const int signs[32] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    double result = signs[mask];
    for (int i = 0; i < 5; i++) if (mask & (1 << i)) result *= k.data[i];
    return result;
}

/* The central unit I has square s = product(K), including s=0. */
static inline vec2
clifford5_center_mul(double s, vec2 a, vec2 b)
{
    return complex_make((const double[]){a.data[0]*b.data[0] + s*a.data[1]*b.data[1],
        a.data[0]*b.data[1] + a.data[1]*b.data[0]});
}

static inline vec2
clifford5_center_square(double s, vec2 a)
{
    return complex_make((const double[]){a.data[0]*a.data[0] + s*a.data[1]*a.data[1], 2*a.data[0]*a.data[1]});
}

static inline complex
clifford5_center(clifford5 a)
{
    return complex_make((const double[]){a.data[0], a.data[31]});
}

static inline clifford5
clifford5_from_center(complex z)
{
    clifford5 result = {{0}};
    result.data[0] = z.data[0];
    result.data[31] = z.data[1];
    return result;
}

static inline clifford5
clifford5_center_scale(vec5 k, clifford5 a, vec2 z)
{
    const double *x = a.data;
    double re = z.data[0], im = z.data[1];
    double k0 = k.data[0];
    double k1 = k.data[1];
    double k01 = k0*k1;
    double k2 = k.data[2];
    double k02 = k0*k2;
    double k12 = k1*k2;
    double k012 = k01*k2;
    double k3 = k.data[3];
    double k03 = k0*k3;
    double k13 = k1*k3;
    double k013 = k01*k3;
    double k23 = k2*k3;
    double k023 = k02*k3;
    double k123 = k12*k3;
    double k0123 = k012*k3;
    double k4 = k.data[4];
    double k04 = k0*k4;
    double k14 = k1*k4;
    double k014 = k01*k4;
    double k24 = k2*k4;
    double k024 = k02*k4;
    double k124 = k12*k4;
    double k0124 = k012*k4;
    double k34 = k3*k4;
    double k034 = k03*k4;
    double k134 = k13*k4;
    double k0134 = k013*k4;
    double k234 = k23*k4;
    double k0234 = k023*k4;
    double k1234 = k123*k4;
    double k01234 = k0123*k4;
    return clifford5_make((const double[]){re*x[0] + k01234*im*x[31], re*x[1] + k1234*im*x[30], re*x[2] - k0234*im*x[29], re*x[3] - k234*im*x[28], re*x[4] + k0134*im*x[27], re*x[5] + k134*im*x[26], re*x[6] - k034*im*x[25], re*x[7] - k34*im*x[24], re*x[8] - k0124*im*x[23], re*x[9] - k124*im*x[22], re*x[10] + k024*im*x[21], re*x[11] + k24*im*x[20], re*x[12] - k014*im*x[19], re*x[13] - k14*im*x[18], re*x[14] + k04*im*x[17], re*x[15] + k4*im*x[16], re*x[16] + k0123*im*x[15], re*x[17] + k123*im*x[14], re*x[18] - k023*im*x[13], re*x[19] - k23*im*x[12], re*x[20] + k013*im*x[11], re*x[21] + k13*im*x[10], re*x[22] - k03*im*x[9], re*x[23] - k3*im*x[8], re*x[24] - k012*im*x[7], re*x[25] - k12*im*x[6], re*x[26] + k02*im*x[5], re*x[27] + k2*im*x[4], re*x[28] - k01*im*x[3], re*x[29] - k1*im*x[2], re*x[30] + k0*im*x[1], re*x[31] + im*x[0]});
}

static inline complex
clifford5_center_product(vec5 k, clifford5 a, clifford5 b)
{
    const double *x = a.data, *y = b.data;
    return complex_make((const double[]){clifford5_scalar_product(k, a, b),
        x[0]*y[31] + x[1]*y[30] - x[2]*y[29] + x[3]*y[28] + x[4]*y[27] - x[5]*y[26] + x[6]*y[25] + x[7]*y[24] - x[8]*y[23] + x[9]*y[22] - x[10]*y[21] - x[11]*y[20] + x[12]*y[19] + x[13]*y[18] - x[14]*y[17] + x[15]*y[16] + x[16]*y[15] - x[17]*y[14] + x[18]*y[13] + x[19]*y[12] - x[20]*y[11] - x[21]*y[10] + x[22]*y[9] - x[23]*y[8] + x[24]*y[7] + x[25]*y[6] - x[26]*y[5] + x[27]*y[4] + x[28]*y[3] - x[29]*y[2] + x[30]*y[1] + x[31]*y[0]});
}

static inline void
clifford5_centered_quartic(vec5 k, clifford5 b2, clifford5 b3, complex p[5])
{
    complex s2 = clifford5_center(b2), s3 = clifford5_center(b3);
    complex s4 = clifford5_center_product(k, b2, b2);
    p[0] = vec2_sub(vec2_scale(clifford5_center_square(clifford5_pseudoscalar_square(k), s2), 2), s4);
    p[1] = vec2_scale(s3, -4./3);
    p[2] = vec2_scale(s2, -2);
    p[3] = complex_from_re(0);
    p[4] = complex_from_re(1);
}

/*
 * Minimal polynomial over the centre span{1, w}, w = e0*e1*e2*e3*e4, which is
 * the clifford1 algebra with parameter w*w = k0*k1*k2*k3*k4. Each element
 * satisfies sum(p[i]*a^i) = 0 for central coefficients p[i], in ascending
 * order, where p[i] stands for p[i][0] + p[i][1]*w and p[4] = 1.
 */
static inline void
clifford5_center_minpoly(vec5 k, clifford5 a, clifford1 p[5])
{
    double s = clifford5_pseudoscalar_square(k);
    complex w = clifford5_center(a), w2 = clifford5_center_square(s, w), q[5];
    clifford5 b = a;
    b.data[0] = b.data[31] = 0;
    clifford5 b2 = clifford5_square(k, b), b3 = clifford5_mul(k, b2, b);
    clifford5_centered_quartic(k, b2, b3, q);
    p[0] = vec2_add(vec2_sub(clifford5_center_square(s, w2), clifford5_center_mul(s, w, q[1])),
        vec2_add(clifford5_center_mul(s, w2, q[2]), q[0]));
    p[1] = vec2_sub(vec2_sub(q[1], clifford5_center_mul(s, vec2_scale(w, 2), q[2])), clifford5_center_mul(s, vec2_scale(w, 4), w2));
    p[2] = vec2_add(vec2_scale(w2, 6), q[2]);
    p[3] = vec2_scale(w, -4);
    p[4] = complex_from_re(1);
}

static inline void
clifford5_minpoly(vec5 k, clifford5 a, double result[9])
{
    double s = clifford5_pseudoscalar_square(k);
    complex p[5];
    clifford5_center_minpoly(k, a, p);
    for (int i = 0; i < 9; i++) result[i] = 0;
    for (int i = 0; i < 5; i++) {
        result[2*i] += p[i].data[0]*p[i].data[0] - s*p[i].data[1]*p[i].data[1];
        for (int j = 0; j < i; j++) {
            result[i + j] += 2*(p[i].data[0]*p[j].data[0] - s*p[i].data[1]*p[j].data[1]);
        }
    }
}

static inline double
clifford5_norm8(vec5 k, clifford5 a)
{
    vec2 p[5];
    clifford5_center_minpoly(k, a, p);
    return p[0].data[0]*p[0].data[0] - clifford5_pseudoscalar_square(k)*p[0].data[1]*p[0].data[1];
}

static inline complex
clifford5_norm(vec5 k, clifford5 a)
{
    return complex_root_real(clifford5_norm8(k, a), 8);
}

static inline clifford5
clifford5_inv(vec5 k, clifford5 a)
{
    vec2 p[5];
    clifford5_center_minpoly(k, a, p);
    double determinant = p[0].data[0]*p[0].data[0] - clifford5_pseudoscalar_square(k)*p[0].data[1]*p[0].data[1];
    vec2 factor = {{-p[0].data[0]/determinant, p[0].data[1]/determinant}};
    clifford5 a2 = clifford5_square(k, a), a3 = clifford5_mul(k, a2, a);
    clifford5 numerator = vec32_add(a3, vec32_add(clifford5_center_scale(k, a2, p[3]),
        vec32_add(clifford5_center_scale(k, a, p[2]), clifford5_from_center(p[1]))));
    return clifford5_center_scale(k, numerator, factor);
}

static inline void
clifford5_eivals(vec5 k, clifford5 a, complex result[8])
{
    double s = clifford5_pseudoscalar_square(k), radius = sqrt(fabs(s));
    vec2 q[5];
    clifford5 b = a;
    b.data[0] = b.data[31] = 0;
    clifford5 b2 = clifford5_square(k, b), b3 = clifford5_mul(k, b2, b);
    clifford5_centered_quartic(k, b2, b3, q);
    if (s < 0) {
        complex p[5], w = {{a.data[0], radius*a.data[31]}};
        for (int i = 0; i < 5; i++) p[i] = complex_make((const double[]){q[i].data[0], radius*q[i].data[1]});
        poly4_complex_coeffs_to_roots(p, result);
        for (int i = 0; i < 4; i++) {
            result[i] = vec2_add(w, result[i]);
            result[i+4] = complex_conjugate(result[i]);
        }
    } else {
        for (int channel = 0; channel < 2; channel++) {
            double r = channel == 0 ? radius : -radius;
            complex p[5];
            for (int i = 0; i < 5; i++) p[i] = complex_from_re(q[i].data[0] + r*q[i].data[1]);
            poly4_complex_coeffs_to_roots(p, result + 4*channel);
            for (int i = 0; i < 4; i++) result[4*channel+i].data[0] += a.data[0] + r*a.data[31];
        }
    }
}

static inline clifford5
clifford5_div(vec5 k, clifford5 a, clifford5 b)
{
    return clifford5_mul(k, a, clifford5_inv(k, b));
}

static inline clifford5
clifford5_real_slice(vec5 k, complex z)
{
    clifford5 result = clifford5_from_re(z.data[0]);
    if (z.data[1] == 0) return result;
    for (int mask = 1; mask < 32; mask++) {
        double square = 1;
        for (int i = 0; i < 5; i++) if (mask & (1 << i)) square *= k.data[i];
        int grade = popcount(mask);
        if ((grade*(grade-1)/2) & 1) square = -square;
        if (square < 0) { result.data[mask] = z.data[1]/sqrt(-square); return result; }
    }
    return clifford5_nan();
}

static inline int
clifford5_unused_generator(clifford5 a)
{
    for (int generator = 0; generator < 5; generator++) {
        int used = 0;
        for (int i = 1; i < 32; i++) {
            if ((i & (1 << generator)) && a.data[i] != 0) { used = 1; break; }
        }
        if (!used) return generator;
    }
    return -1;
}

static inline clifford5
clifford5_apply_subalgebra(vec5 k, clifford5 a, int generator,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    vec4 subk;
    clifford4 suba;
    clifford5 result = {{0}};
    int low = (1 << generator) - 1, j = 0;
    for (int i = 0; i < 5; i++) if (i != generator) subk.data[j++] = k.data[i];
    for (int i = 0; i < 16; i++) suba.data[i] = a.data[(i & low) | ((i & ~low) << 1)];
    clifford4 value = clifford4_apply(subk, suba, f, deriv1, deriv2, deriv3, cf, cderiv1);
    if (isnan(value.data[0])) return clifford5_nan();
    for (int i = 0; i < 16; i++) result.data[(i & low) | ((i & ~low) << 1)] = value.data[i];
    return result;
}

static inline complex
clifford5_complex_center(double radius, clifford5 a)
{
    return complex_make((const double[]){a.data[0], radius*a.data[31]});
}

static inline clifford5
clifford5_from_complex_center(double radius, complex z)
{
    return clifford5_from_center(complex_make((const double[]){z.data[0], z.data[1]/radius}));
}

static inline clifford5
clifford5_complex_scale(vec5 k, double radius, clifford5 a, complex z)
{
    return clifford5_center_scale(k, a, complex_make((const double[]){z.data[0], z.data[1]/radius}));
}

static inline complex
clifford5_complex_coordinate(vec5 k, double radius, clifford5 a, int i)
{
    return complex_make((const double[]){a.data[i], radius*a.data[31 ^ i]/clifford5_pseudoscalar_product(k, i)});
}

static inline void
clifford5_complex_quartic(vec5 k, double radius, clifford5 b2, clifford5 b3, complex p[5])
{
    vec2 q[5];
    clifford5_centered_quartic(k, b2, b3, q);
    for (int i = 0; i < 5; i++) p[i] = complex_make((const double[]){q[i].data[0], radius*q[i].data[1]});
}

static inline int
clifford5_complex_quadratic_relation(vec5 k, double radius, clifford5 b, clifford5 b2, complex *linear)
{
    int pivot = 1;
    double size = 0;
    for (int i = 2; i < 16; i++) {
        if (poly_complex_abs(clifford5_complex_coordinate(k, radius, b, i)) > poly_complex_abs(clifford5_complex_coordinate(k, radius, b, pivot))) pivot = i;
    }
    *linear = complex_div(clifford5_complex_coordinate(k, radius, b2, pivot), clifford5_complex_coordinate(k, radius, b, pivot));
    if (!hermite_complex_is_finite(*linear)) return 0;
    clifford5 product = clifford5_complex_scale(k, radius, b, *linear);
    for (int i = 0; i < 32; i++) size = fmax(size, fabs(b2.data[i]) + fabs(product.data[i]));
    for (int i = 1; i < 31; i++) {
        if (fabs(b2.data[i] - product.data[i]) > 32*DBL_EPSILON*size) return 0;
    }
    return 1;
}

static inline int
clifford5_complex_cubic_relation(vec5 k, double radius, clifford5 b, clifford5 b2, clifford5 b3, complex p[4])
{
    int first = 1, second = 1;
    double largest = 0, bound = 0;
    for (int i = 2; i < 16; i++) {
        if (poly_complex_abs(clifford5_complex_coordinate(k, radius, b, i)) > poly_complex_abs(clifford5_complex_coordinate(k, radius, b, first))) first = i;
    }
    complex ratio = complex_div(clifford5_complex_coordinate(k, radius, b2, first), clifford5_complex_coordinate(k, radius, b, first));
    for (int i = 1; i < 16; i++) {
        double residual = poly_complex_abs(vec2_sub(clifford5_complex_coordinate(k, radius, b2, i), complex_mul(ratio, clifford5_complex_coordinate(k, radius, b, i))));
        if (residual > largest) { largest = residual; second = i; }
    }
    if (largest == 0) return 0;
    complex quadratic = complex_div(vec2_sub(clifford5_complex_coordinate(k, radius, b3, second),
        complex_mul(complex_div(clifford5_complex_coordinate(k, radius, b3, first), clifford5_complex_coordinate(k, radius, b, first)), clifford5_complex_coordinate(k, radius, b, second))),
        vec2_sub(clifford5_complex_coordinate(k, radius, b2, second), complex_mul(ratio, clifford5_complex_coordinate(k, radius, b, second))));
    complex linear = complex_div(vec2_sub(clifford5_complex_coordinate(k, radius, b3, first),
        complex_mul(quadratic, clifford5_complex_coordinate(k, radius, b2, first))), clifford5_complex_coordinate(k, radius, b, first));
    if (!hermite_complex_is_finite(linear) || !hermite_complex_is_finite(quadratic)) return 0;
    clifford5 q = clifford5_complex_scale(k, radius, b2, quadratic), l = clifford5_complex_scale(k, radius, b, linear);
    for (int i = 0; i < 32; i++) bound = fmax(bound, fabs(b3.data[i]) + fabs(q.data[i]) + fabs(l.data[i]));
    for (int i = 1; i < 31; i++) {
        if (fabs(b3.data[i] - q.data[i] - l.data[i]) > 64*DBL_EPSILON*bound) return 0;
    }
    p[0] = vec2_sub(complex_mul(quadratic, clifford5_complex_center(radius, b2)), clifford5_complex_center(radius, b3));
    p[1] = vec2_neg(linear);
    p[2] = vec2_neg(quadratic);
    p[3] = complex_from_re(1);
    return 1;
}

static inline clifford5
clifford5_apply_complex_center(vec5 k, double radius, clifford5 a,
    double (*f)(double), complex (*cf)(complex), complex (*cderiv1)(complex),
    complex (*cderiv2)(complex), complex (*cderiv3)(complex))
{
    complex w = clifford5_complex_center(radius, a), linear, coefficients[4], p[5];
    clifford5 b = a;
    b.data[0] = b.data[31] = 0;
    if (!vec32_any_nonzero(b)) return clifford5_from_complex_center(radius, cf(w));
    clifford5 b2 = clifford5_square(k, b), b3;
    int ok;
    if (clifford5_complex_quadratic_relation(k, radius, b, b2, &linear)) {
        complex shift = vec2_scale(linear, .5);
        complex root = complex_sqrt(vec2_add(clifford5_complex_center(radius, b2), complex_square(shift)));
        complex nodes[2] = {vec2_sub(shift, root), vec2_add(shift, root)};
        ok = hermite_complex_interpolate(2, w, nodes,
            f, cf, cderiv1, cderiv2, cderiv3, coefficients);
        if (!ok) return clifford5_nan();
        return vec32_add(clifford5_from_complex_center(radius, coefficients[0]), clifford5_complex_scale(k, radius, b, coefficients[1]));
    }
    b3 = clifford5_mul(k, b2, b);
    if (clifford5_complex_cubic_relation(k, radius, b, b2, b3, p)) {
        complex roots[3];
        poly3_complex_coeffs_to_roots(p, roots);
        ok = hermite_complex_interpolate(3, w, roots,
            f, cf, cderiv1, cderiv2, cderiv3, coefficients);
    } else {
        clifford5_complex_quartic(k, radius, b2, b3, p);
        ok = poly4_complex_apply(p, w, f, cf, cderiv1, cderiv2, cderiv3, coefficients);
    }
    if (!ok) return clifford5_nan();
    return vec32_add(clifford5_from_complex_center(radius, coefficients[0]),
        vec32_add(clifford5_complex_scale(k, radius, b, coefficients[1]),
            vec32_add(clifford5_complex_scale(k, radius, b2, coefficients[2]), clifford5_complex_scale(k, radius, b3, coefficients[3]))));
}

/* Project I to +/-sqrt(s), leaving the first four generators. */
static inline clifford5
clifford5_apply_split_center(vec5 k, double radius, clifford5 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex))
{
    vec4 subk = {{k.data[0], k.data[1], k.data[2], k.data[3]}};
    clifford4 plus, minus;
    clifford5 result;
    for (int i = 0; i < 16; i++) {
        double offset = radius*a.data[31 ^ i]/clifford5_pseudoscalar_product(k, i);
        plus.data[i] = a.data[i] + offset;
        minus.data[i] = a.data[i] - offset;
    }
    plus = clifford4_apply(subk, plus, f, deriv1, deriv2, deriv3, cf, cderiv1);
    minus = clifford4_apply(subk, minus, f, deriv1, deriv2, deriv3, cf, cderiv1);
    if (isnan(plus.data[0]) || isnan(minus.data[0])) return clifford5_nan();
    for (int i = 0; i < 16; i++) {
        result.data[i] = plus.data[i]/2 + minus.data[i]/2;
        result.data[31 ^ i] = (plus.data[i]/2 - minus.data[i]/2)*clifford5_pseudoscalar_product(k, i)/radius;
    }
    return result;
}

static inline int
clifford5_quadratic_relation(clifford5 b, clifford5 b2, double *linear)
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
clifford5_cubic_relation(clifford5 b, clifford5 b2, clifford5 b3, double p[4])
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

/* In a dual centre the real quartic annihilates the semisimple quotient.
 * Use it directly when it also annihilates this element, avoiding derivatives
 * at a semisimple zero eigenvalue that has no nilpotent part. */
static inline int
clifford5_quotient_relation(vec5 k, clifford5 b, clifford5 b2, const double p[5])
{
    clifford5 b4 = clifford5_square(k, b2);
    double error = 0, bound = fabs(p[0]);
    for (int i = 0; i < 32; i++) {
        double value = b4.data[i] + p[2]*b2.data[i] + p[1]*b.data[i] + (i == 0 ? p[0] : 0);
        error = fmax(error, fabs(value));
        bound = fmax(bound, fabs(b4.data[i]) + fabs(p[2]*b2.data[i]) + fabs(p[1]*b.data[i]));
    }
    return error <= 64*DBL_EPSILON*bound;
}

/* Degree-eight real Hermite calculus for a dual centre or close split channels. */
static inline clifford5
clifford5_apply_real_spectrum(vec5 k, clifford5 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    double (*deriv4)(double), double (*deriv5)(double),
    double (*deriv6)(double), double (*deriv7)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex),
    complex (*cderiv2)(complex), complex (*cderiv3)(complex))
{
    double w = a.data[0], linear, coefficients[8], cubic[4];
    clifford5 b = a;
    b.data[0] = 0;
    clifford5 b2 = clifford5_square(k, b), b3;
    int ok;
    if (clifford5_quadratic_relation(b, b2, &linear)) {
        double shift = linear/2, output[4];
        complex root = complex_sqrt_real(b2.data[0] + shift*shift);
        complex nodes[2] = {complex_ladd(vec2_neg(root), shift), complex_ladd(root, shift)};
        ok = hermite_interpolate(2, w, nodes, f, deriv1, deriv2, deriv3, cf, cderiv1, output);
        if (!ok) return clifford5_nan();
        clifford5 result = vec32_scale(b, output[1]);
        result.data[0] = output[0];
        return result;
    }
    b3 = clifford5_mul(k, b2, b);
    if (clifford5_cubic_relation(b, b2, b3, cubic)) {
        complex nodes[3];
        double output[4];
        poly3_coeffs_to_roots(cubic, nodes);
        ok = hermite_interpolate(3, w, nodes, f, deriv1, deriv2, deriv3, cf, cderiv1, output);
        if (!ok) return clifford5_nan();
        clifford5 result = vec32_add(vec32_scale(b, output[1]), vec32_scale(b2, output[2]));
        result.data[0] += output[0];
        return result;
    }
    complex nodes[8];
    if (clifford5_pseudoscalar_square(k) == 0) {
        double s2 = b2.data[0], s3 = b3.data[0], s4 = clifford5_scalar_product(k, b2, b2);
        double q[5] = {2*s2*s2-s4, -4*s3/3, -2*s2, 0, 1};
        if (clifford5_quotient_relation(k, b, b2, q)) {
            double output[4];
            if (!poly4_apply(q, w, f, deriv1, deriv2, deriv3, cf, cderiv1, output)) return clifford5_nan();
            clifford5 result = vec32_add(vec32_scale(b, output[1]),
                vec32_add(vec32_scale(b2, output[2]), vec32_scale(b3, output[3])));
            result.data[0] += output[0];
            return result;
        }
        complex p[5];
        for (int i = 0; i < 5; i++) p[i] = complex_from_re(q[i]);
        poly4_complex_coeffs_to_roots(p, nodes);
        for (int i = 0; i < 4; i++) nodes[4+i] = nodes[i];
    } else {
        clifford5_eivals(k, b, nodes);
        for (int i = 0; i < 8; i++) {
            if (hermite_function_distance(complex_ladd(nodes[i], w), f) == 0) {
                return clifford5_apply_split_center(k, sqrt(clifford5_pseudoscalar_square(k)),
                    a, f, deriv1, deriv2, deriv3, cf, cderiv1);
            }
        }
    }
    double (*derivatives[7])(double) = {deriv1, deriv2, deriv3, deriv4, deriv5, deriv6, deriv7};
    complex (*cderivatives[3])(complex) = {cderiv1, cderiv2, cderiv3};
    if (!hermite8_interpolate(8, w, nodes, f, derivatives, cf, cderivatives, coefficients)) return clifford5_nan();
    clifford5 result, power = b3;
    for (int i = 0; i < 32; i++) result.data[i] = coefficients[1]*b.data[i] + coefficients[2]*b2.data[i] + coefficients[3]*b3.data[i];
    result.data[0] += coefficients[0];
    for (int degree = 4; degree < 8; degree++) {
        power = clifford5_mul(k, power, b);
        result = vec32_add(result, vec32_scale(power, coefficients[degree]));
    }
    return result;
}

/* The centre selects complex, split or dual calculus. */
static inline clifford5
clifford5_apply(vec5 k, clifford5 a,
    double (*f)(double), double (*deriv1)(double),
    double (*deriv2)(double), double (*deriv3)(double),
    double (*deriv4)(double), double (*deriv5)(double),
    double (*deriv6)(double), double (*deriv7)(double),
    complex (*cf)(complex), complex (*cderiv1)(complex),
    complex (*cderiv2)(complex), complex (*cderiv3)(complex))
{
    if (vec32_is_real(a)) return clifford5_real_slice(k, cf(complex_from_re(a.data[0])));
    int generator = clifford5_unused_generator(a);
    if (generator >= 0) return clifford5_apply_subalgebra(k, a, generator, f, deriv1, deriv2, deriv3, cf, cderiv1);
    double s = clifford5_pseudoscalar_square(k);
    if (s < 0) return clifford5_apply_complex_center(k, sqrt(-s), a, f, cf, cderiv1, cderiv2, cderiv3);
    if (s > 0) {
        double radius = sqrt(s), offset = 0, size = 1;
        for (int i = 0; i < 16; i++) {
            offset = fmax(offset, fabs(radius*a.data[31 ^ i]/clifford5_pseudoscalar_product(k, i)));
            size = fmax(size, fabs(a.data[i]));
        }
        if (offset >= 1e-4*size) return clifford5_apply_split_center(k, radius, a, f, deriv1, deriv2, deriv3, cf, cderiv1);
    }
    return clifford5_apply_real_spectrum(k, a, f, deriv1, deriv2, deriv3, deriv4, deriv5, deriv6, deriv7, cf, cderiv1, cderiv2, cderiv3);
}

static inline clifford5
clifford5_exp(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_exp,
        double_exp_deriv1, double_exp_deriv2, double_exp_deriv3, double_exp_deriv4, double_exp_deriv5, double_exp_deriv6, double_exp_deriv7,
        complex_exp, complex_exp_deriv1, complex_exp_deriv2, complex_exp_deriv3);
}

static inline clifford5
clifford5_expm1(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_expm1,
        double_expm1_deriv1, double_expm1_deriv2, double_expm1_deriv3, double_expm1_deriv4, double_expm1_deriv5, double_expm1_deriv6, double_expm1_deriv7,
        complex_expm1, complex_expm1_deriv1, complex_expm1_deriv2, complex_expm1_deriv3);
}

static inline clifford5
clifford5_log(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_log,
        double_log_deriv1, double_log_deriv2, double_log_deriv3, double_log_deriv4, double_log_deriv5, double_log_deriv6, double_log_deriv7,
        complex_log, complex_log_deriv1, complex_log_deriv2, complex_log_deriv3);
}

static inline clifford5
clifford5_log1p(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_log1p,
        double_log1p_deriv1, double_log1p_deriv2, double_log1p_deriv3, double_log1p_deriv4, double_log1p_deriv5, double_log1p_deriv6, double_log1p_deriv7,
        complex_log1p, complex_log1p_deriv1, complex_log1p_deriv2, complex_log1p_deriv3);
}

static inline clifford5
clifford5_sqrt(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_sqrt,
        double_sqrt_deriv1, double_sqrt_deriv2, double_sqrt_deriv3, double_sqrt_deriv4, double_sqrt_deriv5, double_sqrt_deriv6, double_sqrt_deriv7,
        complex_sqrt, complex_sqrt_deriv1, complex_sqrt_deriv2, complex_sqrt_deriv3);
}

static inline clifford5
clifford5_sin(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_sin,
        double_sin_deriv1, double_sin_deriv2, double_sin_deriv3, double_sin_deriv4, double_sin_deriv5, double_sin_deriv6, double_sin_deriv7,
        complex_sin, complex_sin_deriv1, complex_sin_deriv2, complex_sin_deriv3);
}

static inline clifford5
clifford5_cos(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_cos,
        double_cos_deriv1, double_cos_deriv2, double_cos_deriv3, double_cos_deriv4, double_cos_deriv5, double_cos_deriv6, double_cos_deriv7,
        complex_cos, complex_cos_deriv1, complex_cos_deriv2, complex_cos_deriv3);
}

static inline clifford5
clifford5_tan(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_tan,
        double_tan_deriv1, double_tan_deriv2, double_tan_deriv3, double_tan_deriv4, double_tan_deriv5, double_tan_deriv6, double_tan_deriv7,
        complex_tan, complex_tan_deriv1, complex_tan_deriv2, complex_tan_deriv3);
}

static inline clifford5
clifford5_asin(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_asin,
        double_asin_deriv1, double_asin_deriv2, double_asin_deriv3, double_asin_deriv4, double_asin_deriv5, double_asin_deriv6, double_asin_deriv7,
        complex_asin, complex_asin_deriv1, complex_asin_deriv2, complex_asin_deriv3);
}

static inline clifford5
clifford5_acos(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_acos,
        double_acos_deriv1, double_acos_deriv2, double_acos_deriv3, double_acos_deriv4, double_acos_deriv5, double_acos_deriv6, double_acos_deriv7,
        complex_acos, complex_acos_deriv1, complex_acos_deriv2, complex_acos_deriv3);
}

static inline clifford5
clifford5_atan(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_atan,
        double_atan_deriv1, double_atan_deriv2, double_atan_deriv3, double_atan_deriv4, double_atan_deriv5, double_atan_deriv6, double_atan_deriv7,
        complex_atan, complex_atan_deriv1, complex_atan_deriv2, complex_atan_deriv3);
}

static inline clifford5
clifford5_sinh(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_sinh,
        double_sinh_deriv1, double_sinh_deriv2, double_sinh_deriv3, double_sinh_deriv4, double_sinh_deriv5, double_sinh_deriv6, double_sinh_deriv7,
        complex_sinh, complex_sinh_deriv1, complex_sinh_deriv2, complex_sinh_deriv3);
}

static inline clifford5
clifford5_cosh(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_cosh,
        double_cosh_deriv1, double_cosh_deriv2, double_cosh_deriv3, double_cosh_deriv4, double_cosh_deriv5, double_cosh_deriv6, double_cosh_deriv7,
        complex_cosh, complex_cosh_deriv1, complex_cosh_deriv2, complex_cosh_deriv3);
}

static inline clifford5
clifford5_tanh(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_tanh,
        double_tanh_deriv1, double_tanh_deriv2, double_tanh_deriv3, double_tanh_deriv4, double_tanh_deriv5, double_tanh_deriv6, double_tanh_deriv7,
        complex_tanh, complex_tanh_deriv1, complex_tanh_deriv2, complex_tanh_deriv3);
}

static inline clifford5
clifford5_asinh(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_asinh,
        double_asinh_deriv1, double_asinh_deriv2, double_asinh_deriv3, double_asinh_deriv4, double_asinh_deriv5, double_asinh_deriv6, double_asinh_deriv7,
        complex_asinh, complex_asinh_deriv1, complex_asinh_deriv2, complex_asinh_deriv3);
}

static inline clifford5
clifford5_acosh(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_acosh,
        double_acosh_deriv1, double_acosh_deriv2, double_acosh_deriv3, double_acosh_deriv4, double_acosh_deriv5, double_acosh_deriv6, double_acosh_deriv7,
        complex_acosh, complex_acosh_deriv1, complex_acosh_deriv2, complex_acosh_deriv3);
}

static inline clifford5
clifford5_atanh(vec5 k, clifford5 a)
{
    return clifford5_apply(k, a, double_atanh,
        double_atanh_deriv1, double_atanh_deriv2, double_atanh_deriv3, double_atanh_deriv4, double_atanh_deriv5, double_atanh_deriv6, double_atanh_deriv7,
        complex_atanh, complex_atanh_deriv1, complex_atanh_deriv2, complex_atanh_deriv3);
}

static inline clifford5
clifford5_radd(double a, clifford5 b)
{
    return clifford5_add(clifford5_from_re(a), b);
}

static inline clifford5
clifford5_ladd(clifford5 a, double b)
{
    return clifford5_add(a, clifford5_from_re(b));
}

static inline clifford5
clifford5_rsub(double a, clifford5 b)
{
    return clifford5_sub(clifford5_from_re(a), b);
}

static inline clifford5
clifford5_lsub(clifford5 a, double b)
{
    return clifford5_sub(a, clifford5_from_re(b));
}

static inline clifford5
clifford5_rmul(vec5 k, double a, clifford5 b)
{
    return clifford5_mul(k, clifford5_from_re(a), b);
}

static inline clifford5
clifford5_lmul(vec5 k, clifford5 a, double b)
{
    return clifford5_mul(k, a, clifford5_from_re(b));
}

static inline clifford5
clifford5_rdiv(vec5 k, double a, clifford5 b)
{
    return clifford5_div(k, clifford5_from_re(a), b);
}

static inline clifford5
clifford5_ldiv(vec5 k, clifford5 a, double b)
{
    return clifford5_div(k, a, clifford5_from_re(b));
}

static inline clifford5
clifford5_ipow(vec5 k, clifford5 a, int64_t b)
{
    clifford5 result = clifford5_one();
    uint64_t exponent = (uint64_t)b;
    if (b == 0) {
        return result;
    }
    if (b == 1) {
        return a;
    }
    if (b == 2) {
        return clifford5_mul(k, a, a);
    }
    if (b == -1) {
        return clifford5_inv(k, a);
    }
    if (b < 0) {
        a = clifford5_inv(k, a);
        /* Unsigned negation also handles INT64_MIN. */
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result = clifford5_mul(k, result, a);
        }
        exponent >>= 1;
        if (exponent != 0) {
            a = clifford5_mul(k, a, a);
        }
    }
    return result;
}

static inline clifford5
clifford5_rpow(vec5 k, clifford5 a, clifford5 b)
{
    return clifford5_exp(k, clifford5_mul(k, clifford5_log(k, a), b));
}

static inline clifford5
clifford5_pow(vec5 k, clifford5 a, clifford5 b)
{
    return clifford5_rpow(k, a, b);
}

static inline clifford5
clifford5_lpow(vec5 k, clifford5 a, clifford5 b)
{
    return clifford5_exp(k, clifford5_mul(k, b, clifford5_log(k, a)));
}

static inline clifford5
clifford5_rpows(vec5 k, double a, clifford5 b)
{
    return clifford5_rpow(k, clifford5_from_re(a), b);
}

static inline clifford5
clifford5_lpows(vec5 k, clifford5 a, double b)
{
    return clifford5_rpow(k, a, clifford5_from_re(b));
}

static inline clifford5
clifford5_mix(vec5 k, clifford5 a, clifford5 b, double t)
{
    clifford5 delta = clifford5_mul(k, clifford5_inv(k, a), b);
    return clifford5_mul(k, a, clifford5_lpows(k, delta, t));
}

static inline void
clifford5_lmat(vec5 k, clifford5 a, double result[ALGPY_CLIFFORD5_N*ALGPY_CLIFFORD5_N])
{
    clifford5 basis = {{0}};
    clifford5 product;
    int row, col;
    for (col = 0; col < ALGPY_CLIFFORD5_N; col++) {
        basis.data[col] = 1;
        product = clifford5_mul(k, a, basis);
        basis.data[col] = 0;
        for (row = 0; row < ALGPY_CLIFFORD5_N; row++) {
            result[row * ALGPY_CLIFFORD5_N + col] = product.data[row];
        }
    }
}

static inline void
clifford5_rmat(vec5 k, clifford5 a, double result[ALGPY_CLIFFORD5_N*ALGPY_CLIFFORD5_N])
{
    clifford5 basis = {{0}};
    clifford5 product;
    int row, col;
    for (row = 0; row < ALGPY_CLIFFORD5_N; row++) {
        basis.data[row] = 1;
        product = clifford5_mul(k, basis, a);
        basis.data[row] = 0;
        for (col = 0; col < ALGPY_CLIFFORD5_N; col++) {
            result[row * ALGPY_CLIFFORD5_N + col] = product.data[col];
        }
    }
}

#endif /* ALGPY_CLIFFORD5_H */
