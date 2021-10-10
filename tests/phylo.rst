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

    >>> seq1 = 'actgggct'
    >>> #         || |||
    >>> seq2 = 'cgtgagct'
