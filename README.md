# algpy
# Build

    pip install .
    python -c "import algpy; print(algpy.hello_world())"

For development use `pip install -e .`.

## Quaternions

    import numpy as np
    from algpy.algebras import quaternion
    quaternion(1, 2, 3, 4) * quaternion(5, 6, 7, 8)

Provides add, subtract, multiply, divide, negate, `exp`, `log`, `norm2` and equality, as scalars and numpy arrays. Access the four components in `(w, x, y, z)` order with `q[i]`.

The core arithmetic is in `src/algpy/algebras/quaternion/quaternion.h` (plain C,
no Python/NumPy dependency); `_quaternion.c` holds the Python/NumPy wrapper,
modelled on NumPy's `_rational_tests.c`.
