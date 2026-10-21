#ifndef ALGPY_VEC16_H
#define ALGPY_VEC16_H

#define ALGPY_VEC16_N 16

typedef struct {
    double data[ALGPY_VEC16_N];
} vec16;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec16, vec16, ALGPY_VEC16_N)

#endif /* ALGPY_VEC16_H */
