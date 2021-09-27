========================
Pairwise Alignment Tests
========================

Pairwise alignment algorithms attempt to find an optimal or
near-optimal alignment between two sequences.

These doctests exercise various components related to this
implementation.


Needleman-Wunsch
================

    >>> from coolseq.align.pairwise import nw_align, print_alignment
    >>> result = nw_align('at', 'aagt')
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = nw_align('gattaca', 'gcatgcu')
    >>> print_alignment(result)
    g-attaca
    | ||  |
    gcatg-cu
    >>> result = nw_align('atgc', 'attgagc')
    >>> print_alignment(result)
    at-g--c
    || |  |
    attgagc


Smith-Waterman-Beyer
====================

The WSBScorer encapsulates scoring logic for the Waterman-Smith-Beyer
algorithm.

    >>> from coolseq.align.pairwise import WSBScorer

The WSB scorer sets a number of scoring options by default.

    >>> scorer = WSBScorer()
    >>> scorer.options
    {'match': 0, 'mismatch': 1, 'gap_start': 1, 'gap_extend': 1}
    >>> scorer.options == scorer.default_options
    True

Scoring options can be set when creating a scorer instance.

    >>> scorer = WSBScorer({'gap_start': 5})
    >>> scorer.options
    {'match': 0, 'mismatch': 1, 'gap_start': 5, 'gap_extend': 1}
    >>> scoring_options = {
    ...     'match': 1,
    ...     'mismatch': -1,
    ...     'gap_start': 0,
    ...     'gap_extend': -1,
    ... }


    >>> from coolseq.align.pairwise import (
    ...     initialize_matrix,
    ...     print_matrix,
    ...     print_arrow_matrix,
    ... )
    >>> scores, arrows = initialize_matrix('ag', 'tagt', WSBScorer())
    >>> print_matrix(scores)
    [0, 2, 3, 4, 5]
    [2, 1, 2, 4, 5]
    [3, 3, 2, 2, 4]
    >>> print_arrow_matrix(arrows)
    [∅ ← ← ← ←]
    [↑ ↖ ↖ ← ←]
    [↑ ↑ ↖ ↖ ←]

Put it all together.

    >>> from coolseq.align.pairwise import wsb_align
    >>> result = wsb_align('at', 'aagt')
    >>> print_alignment(result)
    a--t
    |  |
    aagt
    >>> result = wsb_align('ag', 'tagt')
    >>> print_alignment(result)
    -ag-
     ||
    tagt
    >>> result = wsb_align('gattaca', 'gcatgcu')
    >>> print_alignment(result)
    gattaca
    |  | |
    gcatgcu
    >>> result = wsb_align('gattaca', 'gcatgcu', {'mismatch': 2})
    >>> print_alignment(result)
    g-attaca
    | ||  |
    gcatg-cu

..    >>> print_alignment(result)
..    a--t
..    |  |
..    aagt
..    >>> result = wsb_align('gattaca', 'gcatgcu')
..    >>> print_alignment(result)
..    g-attaca
..    | ||  |
..    gcatg-cu
..    >>> result = wsb_align('atgc', 'attgagc')
..    >>> print_alignment(result)
..    at-g--c
..    || |  |
..    attgagc
