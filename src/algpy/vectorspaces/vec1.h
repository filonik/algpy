#ifndef ALGPY_VEC1_H
#define ALGPY_VEC1_H

#define ALGPY_VEC1_N 1

typedef struct {
    double data[ALGPY_VEC1_N];
} vec1;

typedef struct {
    unsigned char data[ALGPY_VEC1_N];
} bvec1;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec1, vec1, ALGPY_VEC1_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec1, vec1, bvec1, ALGPY_VEC1_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec1, bvec1, ALGPY_VEC1_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec1, vec1, bvec1, bvec1)

#endif /* ALGPY_VEC1_H */