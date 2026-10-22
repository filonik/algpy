#ifndef ALGPY_VEC2_H
#define ALGPY_VEC2_H

#define ALGPY_VEC2_N 2

typedef struct {
    double data[ALGPY_VEC2_N];
} vec2;

typedef struct {
    unsigned char data[ALGPY_VEC2_N];
} bvec2;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec2, vec2, ALGPY_VEC2_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec2, vec2, bvec2, ALGPY_VEC2_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec2, bvec2, ALGPY_VEC2_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec2, vec2, bvec2, bvec2)

#endif /* ALGPY_VEC2_H */
