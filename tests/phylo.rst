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

We'll perform WPGMA on the following example data. Tree data is
returned both in a tuple representation as well as a Scipy linkage
representation which can be used for building dendrograms.

In the tuple representation, a branch is represented as a triple of
(HEIGHT, L_TREE, R_TREE), and a leaf is represented as a pair of
(TAXA_INDEX, TAXA_LABEL).

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
    (5.0,
      (3.0,
        (0, 'A'),
        (1, 'B')),
      (3.0,
        (4, 'E'),
        (1.0,
          (2, 'C'),
          (3, 'D'))))
    >>> print(links)
    [[0. 1. 3. 2.]
     [2. 3. 1. 2.]
     [4. 6. 3. 3.]
     [5. 7. 5. 5.]]

Here's another example using the primate data.

    >>> example2 = [
    ...     [0, 79, 92, 144, 162],
    ...     [79, 0, 95, 154, 169],
    ...     [92, 95, 0, 150, 169],
    ...     [144, 154, 150, 0, 169],
    ...     [162, 169, 169, 169, 0],
    ... ]
    >>> names2 = ['Human', 'Chimp', 'Gorilla', 'Orang-utan', 'Gibbon']
    >>> tree, links, names = wpgma(example2, names2)
    >>> tree
    (84.0625,
      (4, 'Gibbon'),
      (74.75,
        (3, 'Orang-utan'),
        (46.75,
          (2, 'Gorilla'),
          (39.5,
            (0, 'Human'),
            (1, 'Chimp')))))
    >>> print(links)
    [[ 0.      1.     39.5     2.    ]
     [ 2.      5.     46.75    3.    ]
     [ 3.      6.     74.75    4.    ]
     [ 4.      7.     84.0625  5.    ]]


UPGMA
=====

Unweighted pair group method with arithmetic mean (UPGMA) is a
clustering method which can be used to relate taxa and build
phylogenetic trees. This method works similarly to WPGMA above, and
returns results in the same formats.

    >>> from coolseq.phylo import upgma
    >>> tree, links, names = upgma(example1, names1)
    >>> tree
    (5.0,
      (3.0,
        (0, 'A'),
        (1, 'B')),
      (3.0,
        (4, 'E'),
        (1.0,
          (2, 'C'),
          (3, 'D'))))
    >>> print(links)
    [[0. 1. 3. 2.]
     [2. 3. 1. 2.]
     [4. 6. 3. 3.]
     [5. 7. 5. 5.]]

Again using the primate data.

    >>> tree, links, names = wpgma(example2, names2)
    >>> tree
    (84.0625,
      (4, 'Gibbon'),
      (74.75,
        (3, 'Orang-utan'),
        (46.75,
          (2, 'Gorilla'),
          (39.5,
            (0, 'Human'),
            (1, 'Chimp')))))
    >>> print(links)
    [[ 0.      1.     39.5     2.    ]
     [ 2.      5.     46.75    3.    ]
     [ 3.      6.     74.75    4.    ]
     [ 4.      7.     84.0625  5.    ]]

Another example which mirrors the working example from the UPGMA
wikipedia page.

    >>> example3 = [
    ...     [0, 17, 21, 31, 23],
    ...     [17, 0, 30, 34, 21],
    ...     [21, 30, 0, 28, 39],
    ...     [31, 34, 28, 0, 43],
    ...     [23, 21, 39, 43, 0],
    ... ]
    >>> tree, links, names = upgma(example3, names1)
    >>> tree
    (16.5,
      (11.0,
        (4, 'E'),
        (8.5,
          (0, 'A'),
          (1, 'B'))),
      (14.0,
        (2, 'C'),
        (3, 'D')))
    >>> print(links)
    [[ 0.   1.   8.5  2. ]
     [ 4.   5.  11.   3. ]
     [ 2.   3.  14.   2. ]
     [ 6.   7.  16.5  5. ]]
