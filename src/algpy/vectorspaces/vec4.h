#ifndef ALGPY_VEC4_H
#define ALGPY_VEC4_H

#define ALGPY_VEC4_N 4

typedef struct {
    double data[ALGPY_VEC4_N];
} vec4;

typedef struct {
    unsigned char data[ALGPY_VEC4_N];
} bvec4;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec4, vec4, ALGPY_VEC4_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec4, vec4, bvec4, ALGPY_VEC4_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec4, bvec4, ALGPY_VEC4_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec4, vec4, bvec4, bvec4)

#endif /* ALGPY_VEC4_H */
