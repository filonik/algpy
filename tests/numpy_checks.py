import numpy as np
from algpy.algebras import quaternion


def _sample(cls):
    # The same leading coefficients for every dimension up to six.
    return cls(*[2.0, 0.1, 0.2, 0.05, 0.03, 0.02][:cls.N])


def check_exp_log_sqrt_power(cls):
    x = _sample(cls)
    assert np.allclose(cls.to_float(np.log(np.exp(x))), cls.to_float(x))
    s = np.sqrt(x)
    assert np.allclose(cls.to_float(s * s), cls.to_float(x))
    assert np.allclose(cls.to_float(np.power(x, 2)), cls.to_float(x * x))
    assert np.allclose(cls.to_float(x.exp()), cls.to_float(np.exp(x)))


def check_quaternion_exp_log():
    q = quaternion(1.0, 0.1, 0.2, 0.3)
    assert np.allclose(quaternion.to_float(np.exp(np.log(q))), quaternion.to_float(q))


def check_pow_operator(cls):
    x = _sample(cls)
    assert np.allclose(cls.to_float(x ** 2), cls.to_float(x * x))
    assert np.allclose(cls.to_float(x ** x), cls.to_float(np.power(x, x)))
    assert np.allclose(cls.to_float(2 ** x), cls.to_float(np.power(cls(2.0), x)))
