"""Phylogeny functions.

"""

from math import log
from numbers import Real
from typing import Any

import numpy as np

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


def is_leaf(tree):
    """Return True if tree is a leaf, False otherwise."""
    return len(tree) == 2


def make_links(tree, links=None):
    """Recursively add links to a tree."""
    count = 0
    index = 0
    if not links:
        links = []
    height, l_tree, r_tree = tree
    if is_leaf(l_tree) and is_leaf(r_tree):
        l_index, _ = l_tree
        r_index, _ = r_tree
        count = 2
        links.insert(0, [l_index, r_index, height, count])
        index = max(l_index, r_index)
    elif is_leaf(l_tree):
        l_index, _ = l_tree
        r_count, r_index, links = make_links(r_tree, links)
        count = r_count + 1
        index = max(l_index, r_index) + 1
        links.append([l_index, index, height, count])
    elif is_leaf(r_tree):
        r_index, _ = r_tree
        l_count, l_index, links = make_links(l_tree, links)
        count = l_count + 1
        links.append([l_index, r_index, height, count])
        index = max(l_index, r_index) + 1
    else:
        l_count, l_index, links = make_links(l_tree, links)
        r_count, r_index, links = make_links(r_tree, links)
        count = l_count + r_count
        index = max(l_index, r_index) + 1
        links.append([index, index + 1, height, count])
        index = index + 1
    return (count, index, links)


def to_linkage(tree):
    _, _, links = make_links(tree)
    return np.array(links)


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
    clusters = [(i, names[i],) for i in range(len(matrix))]
    while len(matrix) > 1:
        matrix, clusters = wpgma_shrink(matrix, clusters)
        # print_matrix(matrix)
        # print(clusters)
    # linkage = make_links(clusters[0])
    root = clusters[0]
    links = to_linkage(root)
    return root, links, names
