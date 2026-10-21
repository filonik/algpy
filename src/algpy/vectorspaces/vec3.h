#ifndef ALGPY_VEC3_H
#define ALGPY_VEC3_H

#define ALGPY_VEC3_N 3

typedef struct {
    double data[ALGPY_VEC3_N];
} vec3;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec3, vec3, ALGPY_VEC3_N)

#endif /* ALGPY_VEC3_H */
