==================
Pairwise Alignment
==================

Pairwise alignment algorithms attempt to find an optimal or
near-optimal alignment between two sequences.


Needleman-Wunsch
================

    >>> from coolseq.align.pairwise import needleman_wunsch, print_alignment
    >>> alignment = needleman_wunsch('at', 'aagt')
    >>> print_alignment(alignment)
    a--t
    |  |
    aagt
    >>> alignment = needleman_wunsch('gattaca', 'gcatgcu')
    >>> print_alignment(alignment)
    g-attaca
    | ||  |
    gcatg-cu
