"""Phylogeny functions.

"""

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


def jc_match(a: str, b: str) -> int:
    """Return 0 if the characters match; 1 otherwise."""
    if a != b:
        result = 1
    else:
        result = 0
    return result


def jc_distance(seq1: str, seq2: str) -> Distance:
    """Return Jukes-Cantor distance between two nucleic acid sequences."""
    alpha = 1
    result = alpha * sum(jc_match(c1, c2) for (c1, c2) in zip(seq1, seq2))
    result += alpha * len(extent(seq1, seq2))
    return result
