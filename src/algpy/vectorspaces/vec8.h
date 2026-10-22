#ifndef ALGPY_VEC8_H
#define ALGPY_VEC8_H

#define ALGPY_VEC8_N 8

typedef struct {
    double data[ALGPY_VEC8_N];
} vec8;

typedef struct {
    unsigned char data[ALGPY_VEC8_N];
} bvec8;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec8, vec8, ALGPY_VEC8_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec8, vec8, bvec8, ALGPY_VEC8_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec8, bvec8, ALGPY_VEC8_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec8, vec8, bvec8, bvec8)

#endif /* ALGPY_VEC8_H */
