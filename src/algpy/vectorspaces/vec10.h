#ifndef ALGPY_VEC10_H
#define ALGPY_VEC10_H

#define ALGPY_VEC10_N 10

typedef struct {
    double data[ALGPY_VEC10_N];
} vec10;

typedef struct {
    unsigned char data[ALGPY_VEC10_N];
} bvec10;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec10, vec10, ALGPY_VEC10_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec10, vec10, bvec10, ALGPY_VEC10_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec10, bvec10, ALGPY_VEC10_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec10, vec10, bvec10, bvec10)

#endif /* ALGPY_VEC10_H */
