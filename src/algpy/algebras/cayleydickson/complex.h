/* Compatibility export; the shared complex arithmetic lives in fields. */
#ifndef ALGPY_CAYLEYDICKSON_COMPLEX_H
#define ALGPY_CAYLEYDICKSON_COMPLEX_H

#include "../../fields/complex.h"

/* Grade projections of the complex algebra, as for complex1. */
static inline complex
complex_grade0(complex a)
{
    const double *x = a.data;
    return complex_make((const double[]){x[0], 0});
}

static inline complex
complex_grade1(complex a)
{
    const double *x = a.data;
    return complex_make((const double[]){0, x[1]});
}

#endif /* ALGPY_CAYLEYDICKSON_COMPLEX_H */
