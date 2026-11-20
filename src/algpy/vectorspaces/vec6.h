#ifndef ALGPY_VEC6_H
#define ALGPY_VEC6_H

#define ALGPY_VEC6_N 6

typedef struct {
    double data[ALGPY_VEC6_N];
} vec6;

typedef struct {
    unsigned char data[ALGPY_VEC6_N];
} bvec6;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec6, vec6, ALGPY_VEC6_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec6, vec6, bvec6, ALGPY_VEC6_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec6, bvec6, ALGPY_VEC6_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec6, vec6, bvec6, bvec6)

#endif /* ALGPY_VEC6_H */
