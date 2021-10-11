"""Phylogeny functions.

"""

from math import log
from numbers import Real
from typing import Any

from coolseq.align.pairwise import print_matrix


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


def generate_names(start=None):
    """Generate names, 'A' to 'Z', for naming clusters."""
    name = start or 'A'
    while ord(name) <= ord('Z'):
        yield name
        name = chr(ord(name) + 1)


def find_closest(matrix: DistanceMatrix) -> tuple[int, int]:
    """return the closest clusters in the distance matrix."""
    closest = (0, 1)
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            k, l = closest
            if matrix[i][j] < matrix[k][l]:
                closest = (i, j)
    return closest


def wpgma_shrink(matrix: DistanceMatrix, clusters: list[tuple[Any]]) -> tuple[DistanceMatrix, dict[int, str]]:
    closest = find_closest(matrix)
    # print(f'Closest: {closest}, {clusters[closest[0]]} {clusters[closest[1]]}')
    k, l = closest
    new_cluster = (matrix[k][l] / 2, clusters[k], clusters[l])
    clusters.pop(l)
    clusters.pop(k)
    clusters.append(new_cluster)
    new_vec = []
    for i in range(len(matrix)):
        if i not in closest:
            distance = (matrix[i][k] + matrix[i][l]) / 2
            new_vec.append(distance)
    new_matrix = []
    for i in range(len(matrix)):
        if i in closest:
            continue
        vec = matrix[i]
        k, l = closest
        vec = vec[:l] + vec[l+1:]
        vec = vec[:k] + vec[k+1:]
        new_matrix.append(vec)
    for i in range(len(new_vec)):
        new_matrix[i].append(new_vec[i])
    new_vec.append(0)
    new_matrix.append(new_vec)
    return (new_matrix, clusters)


def wpgma(matrix: DistanceMatrix, names: list[str]) -> None:
    """Build phylogenetic tree using WPGMA."""
    result = None
    clusters = [(names[i],) for i in range(len(matrix))]
    while len(matrix) > 1:
        matrix, clusters = wpgma_shrink(matrix, clusters)
        # print_matrix(matrix)
        # print(clusters)
    return clusters
