from typing import TypeVar

T = TypeVar("T")

def id_(x: T) -> T:
    """Identity function."""
    return x
