#ifndef ALGPY_COMPLEX1_H
#define ALGPY_COMPLEX1_H

#include "../cayleydickson/complex.h"

typedef complex complex1;

#define ALGPY_COMPLEX1_N ALGPY_COMPLEX_N

#define complex1_make complex_make
#define complex1_from_re complex_from_re
#define complex1_to_re complex_to_re
#define complex1_one complex_one
#define complex1_mul complex_mul
#define complex1_conjugate complex_conjugate
#define complex1_norm2 complex_norm2
#define complex1_inv complex_inv
#define complex1_div complex_div
#define complex1_basis complex_basis
#define complex1_coeffs complex_coeffs
#define complex1_minpoly complex_minpoly
#define complex1_eivals complex_eivals
#define complex1_re complex_re
#define complex1_im complex_im
#define complex1_radd complex_radd
#define complex1_ladd complex_ladd
#define complex1_rsub complex_rsub
#define complex1_lsub complex_lsub
#define complex1_rmul complex_rmul
#define complex1_lmul complex_lmul
#define complex1_rdiv complex_rdiv
#define complex1_ldiv complex_ldiv
#define complex1_exp complex_exp
#define complex1_log complex_log
#define complex1_pow complex_pow
#define complex1_cos complex_cos
#define complex1_sin complex_sin
#define complex1_tan complex_tan
#define complex1_cosh complex_cosh
#define complex1_sinh complex_sinh
#define complex1_tanh complex_tanh
#define complex1_sqrt complex_sqrt
#define complex1_sqrt_real complex_sqrt_real
#define complex1_root_real complex_root_real

#endif /* ALGPY_COMPLEX1_H */