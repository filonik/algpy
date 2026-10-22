#ifndef ALGPY_VEC16_H
#define ALGPY_VEC16_H

#define ALGPY_VEC16_N 16

typedef struct {
    double data[ALGPY_VEC16_N];
} vec16;

typedef struct {
    unsigned char data[ALGPY_VEC16_N];
} bvec16;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec16, vec16, ALGPY_VEC16_N)
ALGPY_DEFINE_VECTOR_COMPARISONS(vec16, vec16, bvec16, ALGPY_VEC16_N)
ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(bvec16, bvec16, ALGPY_VEC16_N)
ALGPY_DEFINE_VECTOR_REDUCTIONS(vec16, vec16, bvec16, bvec16)

#endif /* ALGPY_VEC16_H */
