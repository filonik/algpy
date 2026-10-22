"""Clifford algebras registered as numpy user-defined dtypes."""
from ._clifford1 import clifford1
from ._clifford2 import clifford2
from ._clifford3 import clifford3
from ._clifford4 import clifford4

_CONSTRUCTORS = {
    1: clifford1,
    2: clifford2,
    3: clifford3,
    4: clifford4,
}


def clifford(*args):
    """Dispatch to ``cliffordN`` where N is the number of arguments."""
    constructor = _CONSTRUCTORS.get(len(args))
    if constructor is None:
        supported = ", ".join(str(n) for n in sorted(_CONSTRUCTORS))
        raise TypeError(
            f"clifford expects {supported} argument(s), got {len(args)}"
        )
    return constructor(*args)


__all__ = [
    "clifford",
    "clifford1",
    "clifford2",
    "clifford3",
    "clifford4",
]
