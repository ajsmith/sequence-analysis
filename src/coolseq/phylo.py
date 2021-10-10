"""Phylogeny functions.

"""

from functools import partial
from numbers import Real


# Type definitions
Distance = Real
DistanceMatrix = list[list[Distance]]


def order_by_length(seq1: str, seq2: str) -> tuple[str, str]:
    """Return the given sequences in ascending order by length."""
    if len(seq1) > len(seq2):
        result = (seq2, seq1)
    else:
        result = (seq1, seq2)
    return result


def extent(seq1: str, seq2: str) -> str:
    """Return the extent of the longer sequence."""
    seq1, seq2 = order_by_length(seq1, seq2)
    n = len(seq1)
    return seq2[n:]


def jc_match(distance: Distance, a: str, b: str) -> Distance:
    """Return 0 if the characters match, or the distance otherwise."""
    if a != b:
        result = distance
    else:
        result = 0
    return result


def jc_distance(seq1: str, seq2: str) -> Distance:
    """Return Jukes-Cantor distance between two nucleic acid sequences."""
    mismatch = 1
    jc_func = partial(jc_match, mismatch)
    result: Distance = 0
    result += sum(jc_func(c1, c2) for (c1, c2) in zip(seq1, seq2))
    result += len(extent(seq1, seq2)) * mismatch
    return result
