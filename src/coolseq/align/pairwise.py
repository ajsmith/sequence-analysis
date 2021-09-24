"""Pairwise alignment functions.

"""
from collections.abc import Generator, Iterable, Sequence
from itertools import islice, tee


# Type definitions for score and arrow matrices
Score = int
Arrow = int
Matrix = list[list[int]]
ScoreMatrix = list[list[Score]]
ArrowMatrix = list[list[Arrow]]


# Arrow value constants
S_ARROW = 0 # Stop
D_ARROW = 1  # Diagonal
T_ARROW = 2  # Top
L_ARROW = 3  # Left


def needleman_wunsch(sequence1: str, sequence2: str) -> list[str]:
    """Return the pairwise alignment found using Needleman-Wunsch."""
    match = 1
    gap = -1
    mismatch = -1
    (scores, arrows) = initialize_matrix(sequence1, sequence2, match, mismatch, gap)
    path = list(trace_path(arrows))
    alignment = build_alignment(sequence1, sequence2, path)
    # score = compute_alignment_score(scores, path)
    return alignment


def initialize_matrix(
        sequence1: str, sequence2: str, match: int, mismatch: int, gap: int
    ) -> tuple[ScoreMatrix, ArrowMatrix]:
    """Return the initialized matrix.

    Here's a small example:

        >>> scores, arrows = initialize_matrix('at', 'aagt', 1, -1, -1)
        >>> print_matrix(scores)
        [0, -1, -2, -3, -4]
        [-1, 1,  0, -1, -2]
        [-2, 0,  0, -1,  0]
        >>> print_matrix(arrows)
        [0, 3, 3, 3, 3]
        [2, 1, 3, 3, 3]
        [2, 2, 1, 3, 1]

    Here's a bigger example taken from the Needleman-Wunsch Wikipedia page:

        >>> scores, arrows = initialize_matrix('gattaca', 'gcatgcu', 1, -1, -1)
        >>> print_matrix(scores)
        [0,  -1, -2, -3, -4, -5, -6, -7]
        [-1,  1,  0, -1, -2, -3, -4, -5]
        [-2,  0,  0,  1,  0, -1, -2, -3]
        [-3, -1, -1,  0,  2,  1,  0, -1]
        [-4, -2, -2, -1,  1,  1,  0, -1]
        [-5, -3, -3, -1,  0,  0,  0, -1]
        [-6, -4, -2, -2, -1, -1,  1,  0]
        [-7, -5, -3, -1, -2, -2,  0,  0]
        >>> print_matrix(arrows)
        [0, 3, 3, 3, 3, 3, 3, 3]
        [2, 1, 3, 3, 3, 3, 3, 3]
        [2, 2, 1, 1, 3, 3, 3, 3]
        [2, 2, 2, 2, 1, 3, 3, 3]
        [2, 2, 2, 2, 2, 1, 3, 3]
        [2, 2, 2, 1, 2, 2, 1, 3]
        [2, 2, 1, 2, 2, 2, 1, 3]
        [2, 2, 2, 1, 3, 2, 2, 1]

    """
    n = len(sequence1)
    m = len(sequence2)
    (scores, arrows) = initialize_matrix_top(m, gap)
    for i in range(1, n + 1):
        # The left-most score is simply the score directly above plus
        # the gap penalty.
        scores.append([scores[i-1][0] + gap])
        # The left-most arrow always points up.
        arrows.append([T_ARROW])
        for j in range(1, m + 1):
            cell = {}
            # The diagonal score is the diagonal neighbor plus the
            # match/mismatch score
            top_left = scores[i-1][j-1] + match_score(i, j, sequence1, sequence2, match, mismatch)
            cell[top_left] = D_ARROW
            # The top score is the top neighbor plus the gap penalty.
            top = scores[i-1][j] + gap
            cell[top] = T_ARROW
            # The left score is the left neighbor plus the gap
            # penalty.
            left = scores[i][j-1] + gap
            cell[left] = L_ARROW
            # Final score is the max of the top-left, top, and left
            # values.
            final_score = max(top_left, top, left)
            scores[i].append(final_score)
            # The arrow points in the direction of the neighbor from
            # where the best score came. This implementation doesn't
            # include branches, but could with a little extra effort.
            arrows[i].append(cell[final_score])
    return (scores, arrows)


def initialize_matrix_top(m_len: int, gap: int) -> tuple[ScoreMatrix, ArrowMatrix]:
    """Return a matrix with only its top initialized.

    `m_length` is the sequence length (of the 2nd sequence). `gap` is
    the gap penalty.

    For example, if the sequence length is 4 and the gap penalty is
    -1, then the top of the matrix looks like:

        >>> scores, arrows = initialize_matrix_top(4, -1)
        >>> scores
        [[0, -1, -2, -3, -4]]
        >>> arrows
        [[0, 3, 3, 3, 3]]

    With a gap penalty of -2, instead we get:

        >>> scores, arrows = initialize_matrix_top(4, -2)
        >>> scores
        [[0, -2, -4, -6, -8]]
        >>> arrows
        [[0, 3, 3, 3, 3]]

    """
    # The top-left score is always zero
    scores = [[0]]
    # The top-left arrow is always the stop arrow
    arrows = [[S_ARROW]]
    for j in range(1, m_len + 1):
        # The next score is simply the previous score on the left plus
        # the gap penalty.
        scores[0].append(scores[0][j-1] + gap)
        # The next arrow always points left.
        arrows[0].append(L_ARROW)
    return (scores, arrows)


def match_score(
        i: int, j: int, sequence1: str, sequence2: str, match: int, mismatch: int
    ) -> int:
    """Return the match score for a position.

    """
    if is_match(i, j, sequence1, sequence2):
        return match
    else:
        return mismatch


def is_match(i: int, j: int, sequence1: str, sequence2: str) -> bool:
    return sequence1[i-1] == sequence2[j-1]


def trace_path(matrix: ArrowMatrix) -> Iterable[tuple[int, int, int]]:
    """Trace the path back through the arrow matrix.

    Here's a small example:

        >>> scores, arrows = initialize_matrix('at', 'aagt', 1, -1, -1)
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

    >>> _, arrows = initialize_matrix('at', 'aagt', 1, -1, -1)
    >>> path = list(trace_path(arrows))
    >>> alignment = build_alignment('at', 'aagt', path)
    >>> print_alignment(alignment)
    a--t
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
    for row in matrix:
        print(row)


def print_alignment(alignment: list[str]) -> None:
    aligned1, aligned2 = alignment
    print(aligned1)
    print(aligned2)
