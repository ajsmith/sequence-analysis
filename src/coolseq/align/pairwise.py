"""Pairwise alignment functions.

This module provides implementations of Needleman-Wunsch and
Waterman-Smith-Beyer algorithms. Clients can call the nw_align() and
wsb_align() functions to calculate alignments. Main logic for these
are found in their scorer implementations, NWScorer and WSBScorer. A
number of support functions are also found here.

"""

from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from itertools import islice
from typing import Optional


# Type definitions
Score = int
Arrow = int
Matrix = list[list[int]]
ScoreMatrix = list[list[Score]]
ArrowMatrix = list[list[Arrow]]
ScoringOptions = dict[str, int]


# Arrow value constants
S_ARROW = 0  # Stop
D_ARROW = 1  # Diagonal
T_ARROW = 2  # Top
L_ARROW = 3  # Left

ARROW_CHAR_MAP = {
    S_ARROW: '∅',
    D_ARROW: '↖',
    T_ARROW: '↑',
    L_ARROW: '←',
}


class Scorer(ABC):
    """Abstract scorer"""

    match = 1
    mismatch = -1
    gap_extend = -1

    @abstractmethod
    def get_top_edges(self, width: int) -> tuple[ScoreMatrix, ArrowMatrix]:
        """Fill the top edge of the score and arrow matrices."""

    @abstractmethod
    def fill_next_row(
            self, scores: ScoreMatrix, arrows: ArrowMatrix, i: int, width: int, sequence1: str, sequence2: str
        ) -> None:
        """Generate the next row of scores and arrows."""

    def match_score(self, i: int, j: int, seq1: str, seq2: str) -> int:
        """Return the match score for a position."""
        if is_match(i, j, seq1, seq2):
            score = self.match
        else:
            score = self.mismatch
        return score


class NWScorer(Scorer):
    """Determine similarity scores for a cell in an alignment score matrix.

    """

    def __init__(self, match: int=1, mismatch: int=-1, gap_extend: int=-1) -> None:
        """Instantiate a scorer."""
        self.match = match
        self.mismatch = mismatch
        self.gap_extend = gap_extend

    def get_top_edges(self, width: int) -> tuple[ScoreMatrix, ArrowMatrix]:
        """Fill the top edge of the score and arrow matrices.

        >>> scorer = NWScorer(1, -1, -1)
        >>> scores, arrows = scorer.get_top_edges(5)
        >>> scores
        [[0, -1, -2, -3, -4]]
        >>> arrows
        [[0, 3, 3, 3, 3]]
        >>> scorer = NWScorer(1, -1, -2)
        >>> scores, arrows = scorer.get_top_edges(4)
        >>> scores
        [[0, -2, -4, -6]]
        >>> arrows
        [[0, 3, 3, 3]]

        """
        scores_top = []
        arrows_top = []
        # The top-left score is always zero and the top-left arrow is
        # always the stop arrow.
        score = 0
        scores_top.append(score)
        arrows_top.append(S_ARROW)
        for _ in range(1, width):
            # The next score is simply the previous score on the left
            # plus the gap penalty, and the next arrow always points
            # left.
            score += self.gap_extend
            scores_top.append(score)
            arrows_top.append(L_ARROW)
        return ([scores_top], [arrows_top])

    def fill_next_row(
            self, scores: ScoreMatrix, arrows: ArrowMatrix, i: int, width: int, sequence1: str, sequence2: str
        ) -> None:
        """Generate the next row of scores and arrows."""
        # The left-most score is simply the score directly above plus
        # the gap penalty.
        score = scores[i-1][0] + self.gap_extend
        # The left-most arrow always points up.
        arrow = T_ARROW
        # Append new lists containing the initial values to the scores
        # and arrows matrices.
        scores.append([score])
        arrows.append([arrow])
        for j in range(1, width):
            score, arrow = self.score(scores, sequence1, sequence2, i, j)
            scores[i].append(score)
            arrows[i].append(arrow)

    def score(self, scores: ScoreMatrix, sequence1: str, sequence2: str, i: int, j: int) -> tuple[Score, Arrow]:
        """Calculate score and arrow values for a cell position."""
        cell = {}
        # The diagonal score is the diagonal neighbor plus the
        # match/mismatch score
        top_left = scores[i-1][j-1] + self.match_score(i, j, sequence1, sequence2)
        cell[top_left] = D_ARROW
        # The top score is the top neighbor plus the gap penalty.
        top = scores[i-1][j] + self.gap_extend
        cell[top] = T_ARROW
        # The left score is the left neighbor plus the gap penalty.
        left = scores[i][j-1] + self.gap_extend
        cell[left] = L_ARROW
        # Final score is the max of the top-left, top, and left
        # values.
        final_score = max(top_left, top, left)
        # The arrow points in the direction of the neighbor from where
        # the best score came. This implementation doesn't include
        # branches, but could with a little extra effort.
        arrow = cell[final_score]
        return (final_score, arrow)


class WSBScorer(Scorer):
    """Distance scorer for Smith-Waterman-Beyer alignments.

    """

    def __init__(self, options: Optional[ScoringOptions] = None):
        opts = self.default_options
        if options:
            opts.update(options)
        self.match = opts['match']
        self.mismatch = opts['mismatch']
        self.gap_start = opts['gap_start']
        self.gap_extend = opts['gap_extend']
        self.options = opts


    @property
    def default_options(self) -> dict[str, int]:
        """Return default scoring options."""
        return {
            'match': 0,
            'mismatch': 1,
            'gap_start': 1,
            'gap_extend': 1,
        }


    def get_top_edges(self, width: int) -> tuple[ScoreMatrix, ArrowMatrix]:
        """Fill the top edge of the score and arrow matrices.

        >>> scorer = WSBScorer()
        >>> scorer.gap_extend
        1
        >>> scores, arrows = scorer.get_top_edges(5)
        >>> scores
        [[0, 2, 3, 4, 5]]
        >>> arrows
        [[0, 3, 3, 3, 3]]
        >>> scorer = WSBScorer({'gap_extend': 2})
        >>> scores, arrows = scorer.get_top_edges(4)
        >>> scores
        [[0, 3, 5, 7]]
        >>> arrows
        [[0, 3, 3, 3]]

        """
        scores : ScoreMatrix = [[]]
        arrows : ArrowMatrix = [[]]
        # The top-left score is always zero and the top-left arrow is
        # always the stop arrow.
        score = 0
        scores[0].append(score)
        arrows[0].append(S_ARROW)
        for j in range(1, width):
            score = self.best_gap_left(scores, 0, j)
            scores[0].append(score)
            arrows[0].append(L_ARROW)
        return (scores, arrows)

    def fill_next_row(
            self, scores: ScoreMatrix, arrows: ArrowMatrix, i: int, width: int, sequence1: str, sequence2: str
        ) -> None:
        """Generate the next row of scores and arrows."""
        # The left-most score is simply the score directly above plus
        # the gap penalty.
        score = self.best_gap_top(scores, i, 0)
        # The left-most arrow always points up.
        arrow = T_ARROW
        # Append new lists containing the initial values to the scores
        # and arrows matrices.
        scores.append([score])
        arrows.append([arrow])
        for j in range(1, width):
            score, arrow = self.score(scores, sequence1, sequence2, i, j)
            scores[i].append(score)
            arrows[i].append(arrow)

    def score(self, scores: ScoreMatrix, sequence1: str, sequence2: str, i: int, j: int) -> tuple[Score, Arrow]:
        """Calculate score and arrow values for a cell position."""
        cell = {}
        # The diagonal score is the diagonal neighbor plus the
        # match/mismatch score
        top_left = scores[i-1][j-1] + self.match_score(i, j, sequence1, sequence2)
        cell[top_left] = D_ARROW
        # The top score is the top neighbor plus the gap penalty.
        top = self.best_gap_top(scores, i, j)
        cell[top] = T_ARROW
        # The left score is the left neighbor plus the gap penalty.
        left = self.best_gap_left(scores, i, j)
        cell[left] = L_ARROW
        # Final score is the min of the top-left, top, and left
        # values.
        final_score = min(top_left, top, left)
        # The arrow points in the direction of the neighbor from where
        # the best score came. This implementation doesn't include
        # branches, but could with a little extra effort.
        arrow = cell[final_score]
        return (final_score, arrow)

    def best_gap_top(self, scores: ScoreMatrix, i: int, j: int) -> int:
        """Find the best gap score top of i, j."""
        return min(scores[i-k][j] + self.gap_penalty(k) for k in range(1, i+1))

    def best_gap_left(self, scores: ScoreMatrix, i: int, j: int) -> int:
        """Find the best gap score left of i, j."""
        return min(scores[i][j-k] + self.gap_penalty(k) for k in range(1, j+1))

    def gap_penalty(self, k: int) -> int:
        """Return affine gap penalty."""
        return self.gap_start + k * self.gap_extend


def nw_align(
        sequence1: str, sequence2: str, opts: Optional[ScoringOptions] = None
    ) -> list[str]:
    """Return the pairwise alignment found using Needleman-Wunsch."""
    opts = opts or {}
    scorer = NWScorer(**opts)
    scores, arrows = initialize_matrix(sequence1, sequence2, scorer)
    path = list(trace_path(arrows))
    alignment = build_alignment(sequence1, sequence2, path)
    # score = compute_alignment_score(scores, path)
    return alignment


def wsb_align(seq1: str, seq2: str, opts: Optional[ScoringOptions] = None) -> list[str]:
    """Align two sequences using Smith-Waterman-Beyer."""
    scorer = WSBScorer(opts)
    scores, arrows = initialize_matrix(seq1, seq2, scorer)
    path = list(trace_path(arrows))
    alignment = build_alignment(seq1, seq2, path)
    # score = compute_alignment_score(scores, path)
    return alignment


def initialize_matrix(
        sequence1: str, sequence2: str, scorer: Scorer,
    ) -> tuple[ScoreMatrix, ArrowMatrix]:
    """Return the initialized matrix."""
    n = len(sequence1)
    m = len(sequence2)
    width = m + 1
    scores, arrows = scorer.get_top_edges(width)
    for i in range(n):
        scorer.fill_next_row(scores, arrows, i + 1, width, sequence1, sequence2)
    return (scores, arrows)


def is_match(i: int, j: int, seq1: str, seq2: str) -> bool:
    """True if the sequences match at positions i, j; False otherwise.

    """
    return seq1[i-1] == seq2[j-1]


def trace_path(arrows: ArrowMatrix) -> Iterable[tuple[int, int, int]]:
    """Trace the path back through the arrow matrix.

    Here's a small example:

    >>> scores, arrows = initialize_matrix(
    ...     'at', 'aagt', NWScorer(1, -1, -1))
    >>> path = trace_path(arrows)
    >>> print(list(path))
    [(2, 4, 1), (1, 3, 3), (1, 2, 3), (1, 1, 1), (0, 0, 0)]

    """
    i = len(arrows) - 1
    j = len(arrows[0]) - 1
    arrow = arrows[i][j]
    yield (i, j, arrow)
    while arrow != S_ARROW:
        if arrow == D_ARROW:
            i -= 1
            j -= 1
        elif arrow == T_ARROW:
            i -= 1
        elif arrow == L_ARROW:
            j -= 1
        arrow = arrows[i][j]
        yield (i, j, arrow)


def build_alignment(sequence1: str, sequence2: str, path: Sequence[tuple[int, int, int]]) -> list[str]:
    """Align two sequences from the arrow path.

    >>> _, arrows = initialize_matrix('at', 'aagt', NWScorer(1, -1, -1))
    >>> path = list(trace_path(arrows))
    >>> alignment = build_alignment('at', 'aagt', path)
    >>> print_alignment(alignment)
    a--t
    |  |
    aagt

    """
    aligned1 = ''
    aligned2 = ''
    for (i, j, arrow) in islice(reversed(path), 1, None):
        i -= 1
        j -= 1
        # print(i, j, arrow)
        if arrow == D_ARROW:
            aligned1 = aligned1 + sequence1[i]
            aligned2 = aligned2 + sequence2[j]
        elif arrow == T_ARROW:
            aligned1 = aligned1 + sequence1[i]
            aligned2 = aligned2 + '-'
        elif arrow == L_ARROW:
            aligned1 = aligned1 + '-'
            aligned2 = aligned2 + sequence2[j]
    return [aligned1, aligned2]


def print_matrix(matrix: Matrix) -> None:
    """Print a matrix."""
    for row in matrix:
        print(row)


def print_alignment(alignment: list[str]) -> None:
    """Print an alignment."""
    aligned1, aligned2 = alignment
    max_line_length = 72
    bars = ''
    for (c1, c2) in zip(aligned1, aligned2):
        if c1 == c2:
            bars += '|'
        else:
            bars += ' '
    while aligned1:
        print(aligned1[:max_line_length])
        print(bars[:max_line_length].rstrip())
        print(aligned2[:max_line_length])
        aligned1 = aligned1[max_line_length:]
        bars = bars[max_line_length:]
        aligned2 = aligned2[max_line_length:]
        if aligned1:
            print()


def print_arrow_matrix(matrix: ArrowMatrix) -> None:
    """Pretty print an arrow matrix."""
    for row in matrix:
        arrow_line = ' '.join(ARROW_CHAR_MAP[e] for e in row)
        print('[' + arrow_line + ']')
