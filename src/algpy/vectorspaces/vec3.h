#ifndef ALGPY_VEC3_H
#define ALGPY_VEC3_H

#define ALGPY_VEC3_N 3

typedef struct {
    double data[ALGPY_VEC3_N];
} vec3;

typedef struct {
    unsigned char data[ALGPY_VEC3_N];
} bvec3;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec3, vec3, ALGPY_VEC3_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec3, vec3, bvec3, ALGPY_VEC3_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec3, bvec3, ALGPY_VEC3_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec3, vec3, bvec3, bvec3)

#endif /* ALGPY_VEC3_H */
