"""Cayley-Dickson algebras registered as numpy user-defined dtypes."""
from ._cayleydickson1 import cayleydickson1
from ._cayleydickson2 import cayleydickson2
from ._cayleydickson3 import cayleydickson3
from ._cayleydickson4 import cayleydickson4

_CONSTRUCTORS = {
    1: cayleydickson1,
    2: cayleydickson2,
    3: cayleydickson3,
    4: cayleydickson4,
}


def cayleydickson(*args):
    """Dispatch to ``cayleydicksonN`` where N is the number of arguments."""
    constructor = _CONSTRUCTORS.get(len(args))
    if constructor is None:
        supported = ", ".join(str(n) for n in sorted(_CONSTRUCTORS))
        raise TypeError(
            f"cayleydickson expects {supported} argument(s), got {len(args)}"
        )
    return constructor(*args)


__all__ = [
    "cayleydickson",
    "cayleydickson1",
    "cayleydickson2",
    "cayleydickson3",
    "cayleydickson4",
]
