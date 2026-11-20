#ifndef ALGPY_VEC32_H
#define ALGPY_VEC32_H

#define ALGPY_VEC32_N 32

typedef struct {
    double data[ALGPY_VEC32_N];
} vec32;

typedef struct {
    unsigned char data[ALGPY_VEC32_N];
} bvec32;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec32, vec32, ALGPY_VEC32_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec32, vec32, bvec32, ALGPY_VEC32_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec32, bvec32, ALGPY_VEC32_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec32, vec32, bvec32, bvec32)

#endif /* ALGPY_VEC32_H */
