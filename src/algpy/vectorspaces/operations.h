#ifndef ALGPY_VECTORSPACE_OPERATIONS_H
#define ALGPY_VECTORSPACE_OPERATIONS_H

#define ALGPY_DEFINE_VECTOR_OPERATIONS(type, prefix, n)                      \
    static inline type prefix##_zero(void)                                   \
    {                                                                        \
        type result = {{0}};                                                 \
        return result;                                                       \
    }                                                                        \
                                                                             \
    static inline int prefix##_nonzero(type a)                               \
    {                                                                        \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            if (a.data[i] != 0) {                                            \
                return 1;                                                    \
            }                                                                \
        }                                                                    \
        return 0;                                                            \
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
    }                                                                        \
                                                                             \
    static inline int prefix##_equal(type a, type b)                         \
    {                                                                        \
        int i;                                                               \
        for (i = 0; i < n; i++) {                                            \
            if (a.data[i] != b.data[i]) {                                    \
                return 0;                                                    \
            }                                                                \
        }                                                                    \
        return 1;                                                            \
    }

#endif /* ALGPY_VECTORSPACE_OPERATIONS_H */
