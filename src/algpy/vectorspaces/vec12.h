#ifndef ALGPY_VEC12_H
#define ALGPY_VEC12_H

#define ALGPY_VEC12_N 12

typedef struct {
    double data[ALGPY_VEC12_N];
} vec12;

typedef struct {
    unsigned char data[ALGPY_VEC12_N];
} bvec12;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec12, vec12, ALGPY_VEC12_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec12, vec12, bvec12, ALGPY_VEC12_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec12, bvec12, ALGPY_VEC12_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec12, vec12, bvec12, bvec12)

#endif /* ALGPY_VEC12_H */
