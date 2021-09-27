========================
Pairwise Alignment Tests
========================

Pairwise alignment algorithms attempt to find an optimal or
near-optimal alignment between two sequences.

These doctests exercise various components related to this
implementation.


Needleman-Wunsch
================

Tests for some of the internals.

    >>> from coolseq.align.pairwise import NWScorer
    >>> from coolseq.align.pairwise import initialize_matrix
    >>> from coolseq.align.pairwise import print_matrix
    >>> from coolseq.align.pairwise import print_arrow_matrix

    >>> scores, arrows = initialize_matrix(
    ...    'at', 'aagt', NWScorer(1, -1, -1))
    >>> print_matrix(scores)
    [0, -1, -2, -3, -4]
    [-1, 1,  0, -1, -2]
    [-2, 0,  0, -1,  0]
    >>> print_arrow_matrix(arrows)
    [∅ ← ← ← ←]
    [↑ ↖ ← ← ←]
    [↑ ↑ ↖ ← ↖]

    >>> scores, arrows = initialize_matrix(
    ...     'gattaca', 'gcatgcu', NWScorer(1, -1, -1))
    >>> print_matrix(scores)
    [0,  -1, -2, -3, -4, -5, -6, -7]
    [-1,  1,  0, -1, -2, -3, -4, -5]
    [-2,  0,  0,  1,  0, -1, -2, -3]
    [-3, -1, -1,  0,  2,  1,  0, -1]
    [-4, -2, -2, -1,  1,  1,  0, -1]
    [-5, -3, -3, -1,  0,  0,  0, -1]
    [-6, -4, -2, -2, -1, -1,  1,  0]
    [-7, -5, -3, -1, -2, -2,  0,  0]
    >>> print_arrow_matrix(arrows)
    [∅ ← ← ← ← ← ← ←]
    [↑ ↖ ← ← ← ← ← ←]
    [↑ ↑ ↖ ↖ ← ← ← ←]
    [↑ ↑ ↑ ↑ ↖ ← ← ←]
    [↑ ↑ ↑ ↑ ↑ ↖ ← ←]
    [↑ ↑ ↑ ↖ ↑ ↑ ↖ ←]
    [↑ ↑ ↖ ↑ ↑ ↑ ↖ ←]
    [↑ ↑ ↑ ↖ ← ↑ ↑ ↖]

    >>> scores, arrows = initialize_matrix(
    ...     'atgc', 'attgagc', NWScorer(1, -1, -1))
    >>> print_matrix(scores)
    [0,  -1, -2, -3, -4, -5, -6, -7]
    [-1,  1,  0, -1, -2, -3, -4, -5]
    [-2,  0,  2,  1,  0, -1, -2, -3]
    [-3, -1,  1,  1,  2,  1,  0, -1]
    [-4, -2,  0,  0,  1,  1,  0,  1]
    >>> print_arrow_matrix(arrows)
    [∅ ← ← ← ← ← ← ←]
    [↑ ↖ ← ← ← ← ← ←]
    [↑ ↑ ↖ ← ← ← ← ←]
    [↑ ↑ ↑ ↖ ↖ ← ← ←]
    [↑ ↑ ↑ ↑ ↑ ↖ ← ↖]


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


Waterman-Smith-Beyer
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

The following example shows matrix initialization for a pair of
sequences using the Waterman-Smith-Beyer algorithm.

    >>> scores, arrows = initialize_matrix('at', 'aagt', WSBScorer())
    >>> print_matrix(scores)
    [0, 2, 3, 4, 5]
    [2, 0, 2, 3, 4]
    [3, 2, 1, 3, 3]
    >>> print_arrow_matrix(arrows)
    [∅ ← ← ← ←]
    [↑ ↖ ← ← ←]
    [↑ ↑ ↖ ← ↖]

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
