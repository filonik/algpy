"""Clifford algebras registered as numpy user-defined dtypes."""
from ._clifford1 import clifford1
from ._clifford2 import clifford2

_CONSTRUCTORS = {
    1: clifford1,
    2: clifford2,
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
]
