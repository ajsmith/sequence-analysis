"""Pairwise alignment functions.

"""


def needleman_wunsch(sequence1: str, sequence2: str) -> list[str]:
    """Return the pairwise alignment found using Needleman-Wunsch."""
    match = 1
    gap = -1
    mismatch = -1
    matrix = initialize_matrix(sequence1, sequence2, match, mismatch, gap)
    return []


def initialize_matrix(
        sequence1: str, sequence2: str, match: int, mismatch: int, gap: int
    ) -> list[list[int]]:
    """Return the initialized matrix.

    Here's an example taken from the Needleman-Wunsch Wikipedia page:

        >>> matrix = initialize_matrix('gattaca', 'gcatgcu', 1, -1, -1)
        >>> print_matrix(matrix)
        [0,  -1, -2, -3, -4, -5, -6, -7]
        [-1,  1,  0, -1, -2, -3, -4, -5]
        [-2,  0,  0,  1,  0, -1, -2, -3]
        [-3, -1, -1,  0,  2,  1,  0, -1]
        [-4, -2, -2, -1,  1,  1,  0, -1]
        [-5, -3, -3, -1,  0,  0,  0, -1]
        [-6, -4, -2, -2, -1, -1,  1,  0]
        [-7, -5, -3, -1, -2, -2,  0,  0]

    """
    n = len(sequence1)
    m = len(sequence2)
    result = initialize_matrix_top(m, gap)
    for i in range(1, n + 1):
        result.append([result[i-1][0] + gap])
        for j in range(1, m + 1):
            top_left = result[i-1][j-1] + match_score(i, j, sequence1, sequence2, match, mismatch)
            top = result[i-1][j] + gap
            left = result[i][j-1] + gap
            score = max(top_left, top, left)
            result[i].append(score)
    return result


def initialize_matrix_top(m_len: int, gap: int) -> list[list[int]]:
    """Return a matrix with only its top initialized.

    `m_length` is the sequence length (of the 2nd sequence). `gap` is
    the gap penalty.

    For example, if the sequence length is 4 and the gap penalty is
    -1, then the top of the matrix looks like:

        >>> initialize_matrix_top(4, -1)
        [[0, -1, -2, -3, -4]]

    With a gap penalty of -2, instead we get:

        >>> initialize_matrix_top(4, -2)
        [[0, -2, -4, -6, -8]]

    """
    result = [[0]]
    for j in range(1, m_len + 1):
        result[0].append(result[0][j-1] + gap)
    return result


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


def print_matrix(matrix: list[list[int]]) -> None:
    for row in matrix:
        print(row)
