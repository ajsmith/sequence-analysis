==========================
Phylogenetic Tree Building
==========================

Phylogenetic trees provide a model of relationships between
taxa. Typically this involves modeling evolutionary relationships
between organisms based on sequence data.


Rates of Change in DNA
======================

The Jukes-Cantor model provides a simple measure of distance between
two sequences.

    >>> from coolseq.phylo import jc_distance

Identical sequences should have no distance between them.

    >>> jc_distance('acgt', 'acgt')
    0

Sequences which differ have a positive, non-zero distance between
them.

    >>> distance = jc_distance('tcgt', 'acgt')
    >>> distance > 0
    True

This distance is the same when the argument order is reversed.

    >>> distance == jc_distance('acgt', 'tcgt')
    True

The distance increases with number of differences.

    >>> jc_distance('tcgt', 'acgt') < jc_distance('ttgt', 'acgt')
    True

Sequences must be equal length (implying they've already been
aligned).

    >>> jc_distance('ac', 'a')
    Traceback (most recent call last):
    ValueError: Sequence lengths differ

Sequences must not be empty.

    >>> jc_distance('', '')
    Traceback (most recent call last):
    ValueError: Empty sequences

Another example.

    >>> seq1 = 'actgggct'
    >>> #         || |||
    >>> seq2 = 'cgtgagct'
    >>> distance = jc_distance(seq1, seq2)
    >>> round(distance, 3)
    0.52


WPGMA
=====

Weighted pair group method with arithmetic mean (WPGMA) is a
clustering method which can be used to relate taxa and build
phylogenetic trees.

We'll perform WPGMA on the following example data.

    >>> from coolseq.phylo import wpgma
    >>> example1 = [
    ...     [0, 6, 10, 10, 10],
    ...     [6, 0, 10, 10, 10],
    ...     [10, 10, 0, 2, 6],
    ...     [10, 10, 2, 0, 6],
    ...     [10, 10, 6, 6, 0],
    ... ]
    >>> names1 = list('ABCDE')
    >>> tree, links, names = wpgma(example1, names1)
    >>> tree
    (5.0, (3.0, (0, 'A'), (1, 'B')), (3.0, (4, 'E'), (1.0, (2, 'C'), (3, 'D'))))
    >>> print(links)
    [[2. 3. 1. 2.]
     [0. 1. 3. 2.]
     [4. 5. 3. 3.]
     [6. 7. 5. 5.]]
