#ifndef ALGPY_VEC4_H
#define ALGPY_VEC4_H

#define ALGPY_VEC4_N 4

typedef struct {
    double data[ALGPY_VEC4_N];
} vec4;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec4, vec4, ALGPY_VEC4_N)

#endif /* ALGPY_VEC4_H */
