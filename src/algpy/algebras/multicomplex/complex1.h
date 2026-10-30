#ifndef ALGPY_COMPLEX1_H
#define ALGPY_COMPLEX1_H

#include "../../fields/complex.h"

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
#define complex1_ipow complex_ipow
#define complex1_rpow complex_rpow
#define complex1_pow complex1_rpow
#define complex1_lpow complex_lpow
#define complex1_cos complex_cos
#define complex1_sin complex_sin
#define complex1_tan complex_tan
#define complex1_cosh complex_cosh
#define complex1_sinh complex_sinh
#define complex1_tanh complex_tanh
#define complex1_sqrt complex_sqrt
#define complex1_sqrt_real complex_sqrt_real
#define complex1_root_real complex_root_real

#define complex1_zero complex_zero
#define complex1_is_real complex_is_real
#define complex1_pos complex_pos
#define complex1_neg complex_neg
#define complex1_add complex_add
#define complex1_sub complex_sub
#define complex1_equal complex_equal
#define complex1_not_equal complex_not_equal
#define complex1_square complex_square
#define complex1_norm complex_norm
#define complex1_acos complex_acos
#define complex1_asin complex_asin
#define complex1_atan complex_atan
#define complex1_acosh complex_acosh
#define complex1_asinh complex_asinh
#define complex1_atanh complex_atanh
#define complex1_expm1 complex_expm1
#define complex1_log1p complex_log1p
#define complex1_rpows complex_rpows
#define complex1_lpows complex_lpows
#define complex1_mix complex_mix

#define complex1_lmat complex_lmat
#define complex1_rmat complex_rmat


#endif /* ALGPY_COMPLEX1_H */
