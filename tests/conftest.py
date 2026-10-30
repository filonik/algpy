"""Common fixtures; each test module explicitly names its algebra."""
import pytest
from algebra_reference import algebra_type

@pytest.fixture
def common_algebra(request):
    cls = algebra_type(request.module.ALGEBRA_NAME)
    if hasattr(cls, "K"):
        cls.K = [-1]*int(request.module.ALGEBRA_NAME[-1])
    return cls

@pytest.fixture(scope="module")
def scipy_functions():
    # An unavailable optional dependency skips only the tests that need it.
    from scipy_reference import FUNCTIONS
    return dict(FUNCTIONS)
