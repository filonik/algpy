#ifndef ALGPY_VEC5_H
#define ALGPY_VEC5_H

#define ALGPY_VEC5_N 5

typedef struct {
    double data[ALGPY_VEC5_N];
} vec5;

typedef struct {
    unsigned char data[ALGPY_VEC5_N];
} bvec5;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec5, vec5, ALGPY_VEC5_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec5, vec5, bvec5, ALGPY_VEC5_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec5, bvec5, ALGPY_VEC5_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec5, vec5, bvec5, bvec5)

#endif /* ALGPY_VEC5_H */
