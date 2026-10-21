#ifndef ALGPY_VEC2_H
#define ALGPY_VEC2_H

#define ALGPY_VEC2_N 2

typedef struct {
    double data[ALGPY_VEC2_N];
} vec2;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec2, vec2, ALGPY_VEC2_N)

#endif /* ALGPY_VEC2_H */
