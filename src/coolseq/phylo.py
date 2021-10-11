"""Phylogeny functions.

"""

from math import log
from numbers import Real


# Type definitions
Distance = Real
DistanceMatrix = list[list[Distance]]


def jc_calc(p_distance: float) -> float:
    """Calculate Jukes-Cantor distance given p-distance."""
    # Make this check because sometimes this formula returns zero as
    # -0.0 and it messes with tests
    if p_distance > 0:
        result = (-3 / 4.0) * log(1 - (4 * p_distance / 3.0))
    else:
        result = 0
    return result


def jc_distance(seq1: str, seq2: str) -> Distance:
    """Return Jukes-Cantor distance between two nucleic acid sequences."""
    if len(seq1) != len(seq2):
        raise ValueError('Sequence lengths differ')

    if len(seq1) == 0:
        raise ValueError('Empty sequence')

    n = len(seq1)
    diffs = sum(c1 != c2 for (c1, c2) in zip(seq1, seq2))
    result = jc_calc(diffs / n)
    return result
