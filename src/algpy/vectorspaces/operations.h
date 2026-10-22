#ifndef ALGPY_VECTORSPACE_OPERATIONS_H
#define ALGPY_VECTORSPACE_OPERATIONS_H

#include <math.h>
#include <stdbool.h>

#define ALGPY_DEFINE_VECTOR_OPERATIONS(type, prefix, n)                      \
    static inline type prefix##_zero(void)                                   \
    {                                                                        \
        type result = {{0}};                                                 \
        return result;                                                       \
    }                                                                        \
                                                                             \
    static inline int prefix##_is_real(type a)                               \
    {                                                                        \
        int i;                                                               \
        for (i = 1; i < n; i++) {                                            \
            if (a.data[i] != 0) {                                            \
                return 0;                                                    \
            }                                                                \
        }                                                                    \
        return 1;                                                            \
    }                                                                        \
                                                                             \
    static inline type prefix##_pos(type a)                                  \
    {                                                                        \
        type result = {{0}};                                                 \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = +a.data[i];                                     \
        }                                                                    \
        return result;                                                       \
    }                                                                        \
                                                                             \
    static inline type prefix##_neg(type a)                                  \
    {                                                                        \
        type result = {{0}};                                                 \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = -a.data[i];                                     \
        }                                                                    \
        return result;                                                       \
    }                                                                        \
                                                                             \
    static inline type prefix##_add(type a, type b)                          \
    {                                                                        \
        type result = {{0}};                                                 \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] + b.data[i];                          \
        }                                                                    \
        return result;                                                       \
    }                                                                        \
                                                                             \
    static inline type prefix##_sub(type a, type b)                          \
    {                                                                        \
        type result = {{0}};                                                 \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] - b.data[i];                          \
        }                                                                    \
        return result;                                                       \
    }                                                                        \
                                                                             \
    static inline type prefix##_scale(type a, double scalar)                 \
    {                                                                        \
        type result = {{0}};                                                 \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] * scalar;                             \
        }                                                                    \
        return result;                                                       \
    }                                                                        \
                                                                             \
    static inline type prefix##_div_scalar(type a, double scalar)            \
    {                                                                        \
        type result = {{0}};                                                 \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] / scalar;                             \
        }                                                                    \
        return result;                                                       \
    }

#define ALGPY_DEFINE_VECTOR_COMPARISONS(type, prefix, btype, n)              \
    static inline btype prefix##_eq(type a, type b)                 \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] == b.data[i];                       \
        }                                                                    \
        return result;                                                       \
    } \
\
    static inline btype prefix##_neq(type a, type b)                 \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] != b.data[i];                       \
        }                                                                    \
        return result;                                                       \
    } \
\
    static inline btype prefix##_lt(type a, type b)                 \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] < b.data[i];                       \
        }                                                                    \
        return result;                                                       \
    } \
\
    static inline btype prefix##_le(type a, type b)                 \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] <= b.data[i];                       \
        }                                                                    \
        return result;                                                       \
    } \
\
    static inline btype prefix##_gt(type a, type b)                 \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] > b.data[i];                       \
        }                                                                    \
        return result;                                                       \
    } \
\
    static inline btype prefix##_ge(type a, type b)                 \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = a.data[i] >= b.data[i];                       \
        }                                                                    \
        return result;                                                       \
    } \
                                                                             \
    static inline btype prefix##_isinf(type a)                              \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = isinf(a.data[i]) ? 1 : 0;                              \
        }                                                                    \
        return result;                                                       \
    } \
                                                                             \
    static inline btype prefix##_isnan(type a)                              \
    {                                                                        \
        btype result = {{0}};                                             \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = isnan(a.data[i]) ? 1 : 0;                              \
        }                                                                    \
        return result;                                                       \
    }

#define ALGPY_DEFINE_BOOLEAN_VECTOR_OPERATIONS(btype, prefix, n)             \
    static inline bool prefix##_all(btype a)                                 \
    {                                                                        \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            if (!a.data[i]) {                                                \
                return false;                                                    \
            }                                                                \
        }                                                                    \
        return true;                                                            \
    }                                                                        \
                                                                             \
    static inline bool prefix##_any(btype a)                                 \
    {                                                                        \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            if (a.data[i]) {                                                 \
                return true;                                                 \
            }                                                                \
        }                                                                    \
        return false;                                                        \
    }                                                                        \
                                                                             \
    static inline btype prefix##_not(btype a)                                \
    {                                                                        \
        btype result = {{0}};                                                \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            result.data[i] = !a.data[i];                                     \
        }                                                                    \
        return result;                                                       \
    }

#define ALGPY_DEFINE_VECTOR_REDUCTIONS(type, prefix, btype, bprefix)         \
    static inline btype prefix##_nonzero(type a)                             \
    {                                                                        \
        return prefix##_neq(a, prefix##_zero());                             \
    }                                                                        \
                                                                             \
    static inline bool prefix##_all_eq(type a, type b)                       \
    {                                                                        \
        return bprefix##_all(prefix##_eq(a, b));                             \
    }                                                                        \
                                                                             \
    static inline bool prefix##_any_nonzero(type a)                          \
    {                                                                        \
        return bprefix##_any(prefix##_nonzero(a));                           \
    }

#endif /* ALGPY_VECTORSPACE_OPERATIONS_H */
