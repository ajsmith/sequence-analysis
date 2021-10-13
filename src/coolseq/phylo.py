"""Phylogeny functions.

"""

from math import log
from numbers import Real
from typing import Any

import numpy as np


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
    """Return the closest clusters in the distance matrix."""
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


def make_links(tree, links, n_clusters):
    """Recursively add links to a tree."""
    count = 0
    height, l_tree, r_tree = tree
    if is_leaf(l_tree) and is_leaf(r_tree):
        l_index, _ = l_tree
        r_index, _ = r_tree
        count = 2
        links.append([l_index, r_index, height, count])
    elif is_leaf(l_tree):
        l_index, _ = l_tree
        r_count, r_index, links = make_links(r_tree, links, n_clusters)
        count = r_count + 1
        links.append([l_index, r_index, height, count])
    elif is_leaf(r_tree):
        r_index, _ = r_tree
        l_count, l_index, links = make_links(l_tree, links, n_clusters)
        count = l_count + 1
        links.append([l_index, r_index, height, count])
    else:
        l_count, l_index, links = make_links(l_tree, links, n_clusters)
        r_count, r_index, links = make_links(r_tree, links, n_clusters)
        count = l_count + r_count
        links.append([l_index, r_index, height, count])
    index = len(links) + n_clusters
    return (count, index, links)


def to_linkage(tree, n_clusters):
    """Create linkage data from a tuple representation of a tree."""
    links = []
    _, _, links = make_links(tree, links, n_clusters)
    return np.array(links)


def wpgma_shrink(matrix: DistanceMatrix, clusters: list[tuple[Any]]) -> tuple[DistanceMatrix, dict[int, str]]:
    """Find the next WPGMA cluster and compact the distance matrix."""
    closest = find_closest(matrix)
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


def wpgma(matrix: DistanceMatrix, names: list[str]) -> tuple:
    """Build phylogenetic tree using WPGMA."""
    n_taxa = len(matrix)
    n_clusters = n_taxa - 1
    clusters = [(i, names[i],) for i in range(len(matrix))]
    while len(matrix) > 1:
        matrix, clusters = wpgma_shrink(matrix, clusters)
    root = clusters[0]
    links = to_linkage(root, n_clusters)
    return root, links, names


def upgma_shrink(matrix, clusters, cluster_sizes):
    """Find the next UPGMA cluster and compact the distance matrix."""
    closest = find_closest(matrix)
    k, l = closest
    new_cluster = (matrix[k][l] / 2, clusters[k], clusters[l])
    k_size = cluster_sizes[clusters[k]]
    l_size = cluster_sizes[clusters[l]]
    clusters.pop(l)
    clusters.pop(k)
    clusters.append(new_cluster)
    cluster_sizes[new_cluster] = k_size + l_size
    new_vec = []
    for i in range(len(matrix)):
        if i not in closest:
            distance = (k_size * matrix[i][k] + l_size * matrix[i][l]) / (k_size + l_size)
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


def upgma(matrix: DistanceMatrix, names: list[str]) -> tuple:
    """Build phylogenetic tree using UPGMA."""
    n_taxa = len(matrix)
    n_clusters = n_taxa - 1
    clusters = [(i, names[i]) for i in range(len(matrix))]
    cluster_sizes = dict((cl, 1) for cl in clusters)
    while len(matrix) > 1:
        matrix, clusters = upgma_shrink(matrix, clusters, cluster_sizes)
    root = clusters[0]
    links = to_linkage(root, n_clusters)
    return root, links, names


def nj_find_smallest(matrix):
    """Find the smallest pair values from the distance matrix."""
    result = (0, 1)
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            k, l = result
            if matrix[i][j] < matrix[k][l]:
                result = (i, j)
    return result


def nj_divergence_vector(matrix):
    """Return the total divergence vector of the distance matrix."""
    return [sum(matrix[i]) for i in range(len(matrix))]


def nj_divergence_matrix(matrix, divergence_vector=None):
    """Return the divergence matrix of a distance matrix."""
    result = []
    n = len(matrix)
    r_vec = divergence_vector or nj_divergence_vector(matrix)
    for i in range(n):
        result.append([])
        for j in range(n):
            if i == j:
                result[i].append(0)
            else:
                result[i].append((n - 2) * matrix[i][j] - r_vec[i] - r_vec[j])
    return result


def nj_shrink(matrix, clusters):
    """Find the next Neighbor Joining cluster and shrink the matrix."""
    n = len(matrix)
    div_vec = nj_divergence_vector(matrix)
    div_matrix = nj_divergence_matrix(matrix, div_vec)
    smallest = nj_find_smallest(div_matrix)
    k, l = smallest
    kl_dist = matrix[k][l]
    k_dist = 0.5 * (kl_dist + (div_vec[k] - div_vec[l]) / (n - 2))
    l_dist = kl_dist - k_dist
    new_cluster = ((k_dist, l_dist), clusters[k], clusters[l])
    clusters.pop(l)
    clusters.pop(k)
    clusters.append(new_cluster)
    new_vec = []
    for i in range(len(matrix)):
        if i not in smallest:
            distance = (matrix[i][k] + matrix[i][l] - kl_dist) / 2
            new_vec.append(distance)
    new_matrix = []
    for i in range(len(matrix)):
        if i in smallest:
            continue
        vec = matrix[i]
        vec = vec[:l] + vec[l+1:]
        vec = vec[:k] + vec[k+1:]
        new_matrix.append(vec)
    for i in range(len(new_vec)):
        new_matrix[i].append(new_vec[i])
    new_vec.append(0)
    new_matrix.append(new_vec)
    return (new_matrix, clusters)


def neighbor_joining(matrix: DistanceMatrix, names: list[str]) -> tuple:
    """Build phylogenetic tree using UPGMA."""
    n_taxa = len(matrix)
    n_clusters = n_taxa - 1
    clusters = [(i, names[i]) for i in range(len(matrix))]
    while len(matrix) > 2:
        matrix, clusters = nj_shrink(matrix, clusters)
    dist = matrix[0][1] / 2
    l_tree, r_tree = clusters
    root = ((dist, dist), l_tree, r_tree)
    # links = to_linkage(root, n_clusters)
    # return root, links, names
    return (root, names)
