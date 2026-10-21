#ifndef ALGPY_VEC8_H
#define ALGPY_VEC8_H

#define ALGPY_VEC8_N 8

typedef struct {
    double data[ALGPY_VEC8_N];
} vec8;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec8, vec8, ALGPY_VEC8_N)

#endif /* ALGPY_VEC8_H */
