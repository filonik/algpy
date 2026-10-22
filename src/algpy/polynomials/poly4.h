#ifndef ALGPY_POLY4_H
#define ALGPY_POLY4_H

#include <math.h>

#include "../algebras/cayleydickson/complex.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Roots of the polynomial by Durand-Kerner iteration; p is lowest degree first, monic. */
static inline void
poly4_solve(const double p[5], complex result[4])
{
    complex z[4];
    int i, j, it;
    for (i = 0; i < 4; i++) {
        double r = 1 + fabs(p[0]) + fabs(p[1]) + fabs(p[2]) + fabs(p[3]);
        double t = 0.4 + 2.0 * M_PI * i / 4;
        z[i] = complex_make(r * cos(t) * 0.5, r * sin(t) * 0.5);
    }
    for (it = 0; it < 500; it++) {
        double delta = 0, scale = 1;
        for (i = 0; i < 4; i++) {
            complex f = complex_from_re(1), den = complex_from_re(1), d;
            for (j = 3; j >= 0; j--) {
                f = complex_mul(f, z[i]);
                f = complex_make(f.data[0] + p[j], f.data[1]);
            }
            for (j = 0; j < 4; j++) {
                if (j != i) {
                    den = complex_mul(den, complex_make(
                        z[i].data[0] - z[j].data[0],
                        z[i].data[1] - z[j].data[1]));
                }
            }
            if (complex_norm2(den) == 0) {
                den = complex_make(1e-12, 0);
            }
            d = complex_div(f, den);
            z[i] = complex_make(z[i].data[0] - d.data[0],
                                 z[i].data[1] - d.data[1]);
            delta += fabs(d.data[0]) + fabs(d.data[1]);
            scale += fabs(z[i].data[0]) + fabs(z[i].data[1]);
        }
        if (delta <= 1e-15 * scale) {
            break;
        }
    }
    for (i = 0; i < 4; i++) {
        result[i] = z[i];
    }
}

#endif /* ALGPY_POLY4_H */
