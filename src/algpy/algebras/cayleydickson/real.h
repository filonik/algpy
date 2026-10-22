/*
 * Shared arithmetic of complex numbers (x + y*i), stored as a vec2.
 *
 * Pure C with no dependency on Python or NumPy. The parameter k is passed
 * explicitly to every operation that depends on it.
 */
#ifndef ALGPY_REAL_H
#define ALGPY_REAL_H

#include <math.h>

#include "../../vectorspaces/vec1.h"

#define ALGPY_REAL_N ALGPY_VEC1_N

typedef vec1 real;

static inline real
real_make(double x)
{
    real r = {{x}};
    return r;
}

static inline real
real_one(void)
{
    return real_make(1);
}

static inline real
real_mul(real x, real y)
{
    return real_make(x.data[0]*y.data[0]);
}

static inline double
real_norm2(real a)
{
    return a.data[0] * a.data[0];
}

static inline double
real_norm(real a)
{
    return sqrt(real_norm2(a));
}

static inline real
real_inv(real a)
{
    double n = real_norm2(a);
    return real_make(a.data[0] / n);
}

static inline real
real_div(real a, real b)
{
    return real_mul(a, real_inv(b));
}

#endif /* ALGPY_REAL_H */