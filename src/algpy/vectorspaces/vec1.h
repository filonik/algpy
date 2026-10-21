#ifndef ALGPY_VEC1_H
#define ALGPY_VEC1_H

#define ALGPY_VEC1_N 1

typedef struct {
    double data[ALGPY_VEC1_N];
} vec1;

#include "operations.h"

ALGPY_DEFINE_VECTOR_OPERATIONS(vec1, vec1, ALGPY_VEC1_N)

#endif /* ALGPY_VEC1_H */