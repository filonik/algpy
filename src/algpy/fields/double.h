#ifndef ALGPY_DOUBLE_H
#define ALGPY_DOUBLE_H

#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define double_expm1 expm1
#define double_exp exp
#define double_log1p log1p
#define double_log log
#define double_pow pow
#define double_cos cos
#define double_sin sin
#define double_tan tan
#define double_acos acos
#define double_asin asin
#define double_atan atan
#define double_cosh cosh
#define double_sinh sinh
#define double_tanh tanh
#define double_acosh acosh
#define double_asinh asinh
#define double_atanh atanh
#define double_sqrt sqrt

static inline double
double_inv(double x)
{
    return 1/x;
}

static inline double
double_square(double x)
{
    return x*x;
}

/* Integer powers by squaring; unsigned magnitude also handles INT64_MIN. */
static inline double
double_ipow(double a, int64_t b)
{
    double result = 1;
    uint64_t exponent = (uint64_t)b;
    if (b < 0) {
        a = 1/a;
        exponent = 0 - exponent;
    }
    while (exponent != 0) {
        if (exponent & 1) {
            result *= a;
        }
        exponent >>= 1;
        if (exponent != 0) {
            a *= a;
        }
    }
    return result;
}

/* Ordinary derivatives f^(N), not Taylor coefficients divided by N!.
 * Evaluate on the smooth real domain of f. Singular points follow IEEE
 * floating-point propagation. */
static inline double
double_exp_deriv1(double x)
{
    return double_exp(x);
}

static inline double
double_exp_deriv2(double x)
{
    return exp(x);
}

static inline double
double_exp_deriv3(double x)
{
    return exp(x);
}

static inline double
double_exp_deriv4(double x)
{
    return exp(x);
}

static inline double
double_expm1_deriv1(double x)
{
    return exp(x);
}

static inline double
double_expm1_deriv2(double x)
{
    return exp(x);
}

static inline double
double_expm1_deriv3(double x)
{
    return exp(x);
}

static inline double
double_expm1_deriv4(double x)
{
    return exp(x);
}

static inline double
double_log_deriv1(double x)
{
    return (1/x);
}

static inline double
double_log_deriv2(double x)
{
    return (-1/x)/x;
}

static inline double
double_log_deriv3(double x)
{
    return ((2/x)/x)/x;
}

static inline double
double_log_deriv4(double x)
{
    return (((-6/x)/x)/x)/x;
}

static inline double
double_log1p_deriv1(double x)
{
    return double_log_deriv1(1+x);
}

static inline double
double_log1p_deriv2(double x)
{
    return double_log_deriv2(1+x);
}

static inline double
double_log1p_deriv3(double x)
{
    return double_log_deriv3(1+x);
}

static inline double
double_log1p_deriv4(double x)
{
    return double_log_deriv4(1+x);
}

static inline double
double_inv_deriv1(double x)
{
    return (-1/x)/x;
}

static inline double
double_inv_deriv2(double x)
{
    return ((2/x)/x)/x;
}

static inline double
double_inv_deriv3(double x)
{
    return (((-6/x)/x)/x)/x;
}

static inline double
double_inv_deriv4(double x)
{
    return ((((24/x)/x)/x)/x)/x;
}

static inline double
double_square_deriv1(double x)
{
    return 2*x;
}

static inline double
double_square_deriv2(double x)
{
    (void)x;
    return 2;
}

static inline double
double_square_deriv3(double x)
{
    (void)x;
    return 0;
}

static inline double
double_square_deriv4(double x)
{
    (void)x;
    return 0;
}

static inline double
double_sqrt_deriv1(double x)
{
    return (0.5/sqrt(x));
}

static inline double
double_sqrt_deriv2(double x)
{
    return (-0.25/x)/sqrt(x);
}

static inline double
double_sqrt_deriv3(double x)
{
    return ((0.375/x)/x)/sqrt(x);
}

static inline double
double_sqrt_deriv4(double x)
{
    return (((-0.9375/x)/x)/x)/sqrt(x);
}

static inline double
double_sin_deriv1(double x)
{
    return (cos(x));
}

static inline double
double_sin_deriv2(double x)
{
    return -sin(x);
}

static inline double
double_sin_deriv3(double x)
{
    return -cos(x);
}

static inline double
double_sin_deriv4(double x)
{
    return sin(x);
}

static inline double
double_cos_deriv1(double x)
{
    return (-sin(x));
}

static inline double
double_cos_deriv2(double x)
{
    return -cos(x);
}

static inline double
double_cos_deriv3(double x)
{
    return sin(x);
}

static inline double
double_cos_deriv4(double x)
{
    return cos(x);
}

static inline double
double_tan_deriv1(double x)
{
    double c = cos(x);
    return (1/c)/c;
}

static inline double
double_tan_deriv2(double x)
{
    return 2*tan(x)*double_tan_deriv1(x);
}

static inline double
double_tan_deriv3(double x)
{
    double t = tan(x);
    return 2*(1+3*t*t)*double_tan_deriv1(x);
}

static inline double
double_tan_deriv4(double x)
{
    double t = tan(x);
    return 8*t*(2+3*t*t)*double_tan_deriv1(x);
}

static inline double
double_asin_deriv1(double x)
{
    return (1/(sqrt(1-x)*sqrt(1+x)));
}

static inline double
double_asin_deriv2(double x)
{
    double d = double_asin_deriv1(x);
    return ((x*d)/(1-x))/(1+x);
}

static inline double
double_asin_deriv3(double x)
{
    double r = (1/(1-x))/(1+x);
    double q = (x/(1-x))/(1+x);
    return (r*r+2*q*q)*double_asin_deriv1(x);
}

static inline double
double_asin_deriv4(double x)
{
    double r = (1/(1-x))/(1+x);
    double q = (x/(1-x))/(1+x);
    return 3*q*(3*r*r+2*q*q)*double_asin_deriv1(x);
}

static inline double
double_acos_deriv1(double x)
{
    return (-1/(sqrt(1-x)*sqrt(1+x)));
}

static inline double
double_acos_deriv2(double x)
{
    double d = double_acos_deriv1(x);
    return ((x*d)/(1-x))/(1+x);
}

static inline double
double_acos_deriv3(double x)
{
    return -double_asin_deriv3(x);
}

static inline double
double_acos_deriv4(double x)
{
    return -double_asin_deriv4(x);
}

static inline double
double_atan_deriv1(double x)
{
    return ((1/hypot(1,x))/hypot(1,x));
}

static inline double
double_atan_deriv2(double x)
{
    double h = hypot(1,x);
    return ((-2*(x/h)/h)/h)/h;
}

static inline double
double_atan_deriv3(double x)
{
    double h = hypot(1, x);
    double r = 1/h, t = x/h;
    return (((2*(3*t*t-r*r)/h)/h)/h)/h;
}

static inline double
double_atan_deriv4(double x)
{
    double h = hypot(1, x);
    double r = 1/h, t = x/h;
    return ((((24*t*(r*r-t*t)/h)/h)/h)/h)/h;
}

static inline double
double_sinh_deriv1(double x)
{
    return cosh(x);
}

static inline double
double_sinh_deriv2(double x)
{
    return sinh(x);
}

static inline double
double_sinh_deriv3(double x)
{
    return cosh(x);
}

static inline double
double_sinh_deriv4(double x)
{
    return sinh(x);
}

static inline double
double_cosh_deriv1(double x)
{
    return sinh(x);
}

static inline double
double_cosh_deriv2(double x)
{
    return cosh(x);
}

static inline double
double_cosh_deriv3(double x)
{
    return sinh(x);
}

static inline double
double_cosh_deriv4(double x)
{
    return cosh(x);
}

static inline double
double_tanh_deriv1(double x)
{
    double t = exp(-fabs(x));
    double s = 2*t/(1+t*t);
    return s*s;
}

static inline double
double_tanh_deriv2(double x)
{
    return -2*tanh(x)*double_tanh_deriv1(x);
}

static inline double
double_tanh_deriv3(double x)
{
    double t = tanh(x);
    return 2*(3*t*t-1)*double_tanh_deriv1(x);
}

static inline double
double_tanh_deriv4(double x)
{
    double t = tanh(x);
    return 8*t*(2-3*t*t)*double_tanh_deriv1(x);
}

static inline double
double_asinh_deriv1(double x)
{
    return (1/hypot(1,x));
}

static inline double
double_asinh_deriv2(double x)
{
    double h = hypot(1,x);
    return (-x/h)*(1/h)/h;
}

static inline double
double_asinh_deriv3(double x)
{
    double h = hypot(1, x);
    double r = 1/h, t = x/h;
    return (((2*t*t-r*r)/h)/h)/h;
}

static inline double
double_asinh_deriv4(double x)
{
    double h = hypot(1, x);
    double r = 1/h, t = x/h;
    return (((3*t*(3*r*r-2*t*t)/h)/h)/h)/h;
}

static inline double
double_acosh_deriv1(double x)
{
    return ((1/sqrt(x-1))/sqrt(x+1));
}

static inline double
double_acosh_deriv2(double x)
{
    return ((-x/(x+1))/(x-1))*double_acosh_deriv1(x);
}

static inline double
double_acosh_deriv3(double x)
{
    double r = (1/(x-1))/(x+1);
    double q = (x/(x-1))/(x+1);
    return (2*q*q+r*r)*double_acosh_deriv1(x);
}

static inline double
double_acosh_deriv4(double x)
{
    double r = (1/(x-1))/(x+1);
    double q = (x/(x-1))/(x+1);
    return -3*q*(2*q*q+3*r*r)*double_acosh_deriv1(x);
}

static inline double
double_atanh_deriv1(double x)
{
    return ((1/(1-x))/(1+x));
}

static inline double
double_atanh_deriv2(double x)
{
    double d = double_atanh_deriv1(x);
    return (2*(x*d))*d;
}

static inline double
double_atanh_deriv3(double x)
{
    double r = double_atanh_deriv1(x);
    double q = (x/(1-x))/(1+x);
    return 2*(r*r+3*q*q)*r;
}

static inline double
double_atanh_deriv4(double x)
{
    double r = double_atanh_deriv1(x);
    double q = (x/(1-x))/(1+x);
    return 24*q*(r*r+q*q)*r;
}

#endif /* ALGPY_DOUBLE_H */
