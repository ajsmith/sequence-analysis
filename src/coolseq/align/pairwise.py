"""Pairwise alignment functions.

"""
from collections.abc import Iterable, Sequence
from itertools import islice


# Type definitions for score and arrow matrices
Score = int
Arrow = int
Matrix = list[list[int]]
ScoreMatrix = list[list[Score]]
ArrowMatrix = list[list[Arrow]]


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


class SimilarityScorer:
    """Determine similarity scores for a cell in an alignment score matrix.

    """

    def __init__(self, match: int=1, mismatch: int=-1, gap_penalty: int=-1) -> None:
        """Instantiate a scorer."""
        self.match = match
        self.mismatch = mismatch
        self.gap_penalty = gap_penalty

    def get_top_edges(self, width: int) -> tuple[ScoreMatrix, ArrowMatrix]:
        """Fill the top edge of the score and arrow matrices.

        >>> scorer = SimilarityScorer(1, -1, -1)
        >>> scores, arrows = scorer.get_top_edges(5)
        >>> scores
        [[0, -1, -2, -3, -4]]
        >>> arrows
        [[0, 3, 3, 3, 3]]
        >>> scorer = SimilarityScorer(1, -1, -2)
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
            score += self.gap_penalty
            scores_top.append(score)
            arrows_top.append(L_ARROW)
        return ([scores_top], [arrows_top])

    def fill_next_row(
            self, scores: ScoreMatrix, arrows: ArrowMatrix, i: int, width: int, sequence1: str, sequence2: str
        ) -> None:
        """Generate the next row of scores and arrows."""
        # The left-most score is simply the score directly above plus
        # the gap penalty.
        score = scores[i-1][0] + self.gap_penalty
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
        top = scores[i-1][j] + self.gap_penalty
        cell[top] = T_ARROW
        # The left score is the left neighbor plus the gap penalty.
        left = scores[i][j-1] + self.gap_penalty
        cell[left] = L_ARROW
        # Final score is the max of the top-left, top, and left
        # values.
        final_score = max(top_left, top, left)
        # The arrow points in the direction of the neighbor from where
        # the best score came. This implementation doesn't include
        # branches, but could with a little extra effort.
        arrow = cell[final_score]
        return (final_score, arrow)

    def match_score(self, i: int, j: int, seq1: str, seq2: str) -> int:
        """Return the match score for a position."""
        if is_match(i, j, seq1, seq2):
            score = self.match
        else:
            score = self.mismatch
        return score


def needleman_wunsch(
        sequence1: str, sequence2: str, match: int=1, mismatch: int=-1, gap: int=-1,
    ) -> list[str]:
    """Return the pairwise alignment found using Needleman-Wunsch."""
    scorer = SimilarityScorer(match, mismatch, gap)
    (scores, arrows) = initialize_matrix(sequence1, sequence2, scorer)
    path = list(trace_path(arrows))
    alignment = build_alignment(sequence1, sequence2, path)
    # score = compute_alignment_score(scores, path)
    return alignment


def initialize_matrix(
        sequence1: str, sequence2: str, scorer: SimilarityScorer,
    ) -> tuple[ScoreMatrix, ArrowMatrix]:
    """Return the initialized matrix.

    Here's a small example:

    >>> scores, arrows = initialize_matrix(
    ...    'at', 'aagt', SimilarityScorer(1, -1, -1))
    >>> print_matrix(scores)
    [0, -1, -2, -3, -4]
    [-1, 1,  0, -1, -2]
    [-2, 0,  0, -1,  0]
    >>> print_arrow_matrix(arrows)
    [∅ ← ← ← ←]
    [↑ ↖ ← ← ←]
    [↑ ↑ ↖ ← ↖]

    Here's a bigger example taken from the Needleman-Wunsch Wikipedia page:

    >>> scores, arrows = initialize_matrix(
    ...     'gattaca', 'gcatgcu', SimilarityScorer(1, -1, -1))
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

    """
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


def trace_path(matrix: ArrowMatrix) -> Iterable[tuple[int, int, int]]:
    """Trace the path back through the arrow matrix.

    Here's a small example:

    >>> scores, arrows = initialize_matrix(
    ...     'at', 'aagt', SimilarityScorer(1, -1, -1))
    >>> path = trace_path(arrows)
    >>> print(list(path))
    [(2, 4, 1), (1, 3, 3), (1, 2, 3), (1, 1, 1), (0, 0, 0)]

    """
    i = len(matrix) - 1
    j = len(matrix[0]) - 1
    arrow = matrix[i][j]
    yield (i, j, arrow)
    while arrow != S_ARROW:
        if arrow == D_ARROW:
            i -= 1
            j -= 1
        elif arrow == T_ARROW:
            i -= 1
        elif arrow == L_ARROW:
            j -= 1
        arrow = matrix[i][j]
        yield (i, j, arrow)


def build_alignment(sequence1: str, sequence2: str, path: Sequence[tuple[int, int, int]]) -> list[str]:
    """Align two sequences from the arrow path.

    >>> _, arrows = initialize_matrix('at', 'aagt', SimilarityScorer(1, -1, -1))
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
    bars = ''
    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i]:
            bars += '|'
        else:
            bars += ' '
    bars = bars.rstrip()
    print(aligned1)
    print(bars)
    print(aligned2)


def print_arrow_matrix(matrix: ArrowMatrix) -> None:
    """Pretty print an arrow matrix."""
    for row in matrix:
        arrow_line = ' '.join(ARROW_CHAR_MAP[e] for e in row)
        print('[' + arrow_line + ']')
