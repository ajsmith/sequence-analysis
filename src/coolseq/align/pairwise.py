"""
Pairwise alignment functions.

"""

def needleman_wunsch(sequence1: str, sequence2: str) -> list[str]:
    """Return the pairwise alignment found using Needleman-Wunsch."""
    match = 1
    gap = -1
    mismatch = -1
    n = len(sequence1)
    m = len(sequence2)
    matrix = initialize_matrix(n, m, match, mismatch, gap)
    print(matrix)
    return []


def initialize_matrix(
        n_len: int, m_len: int, match: int, mismatch: int, gap: int
    ) -> list[list[int]]:
    """Return the initialized matrix.

    """
    result = initialize_matrix_top(m_len, gap)
    for i in range(1, n_len + 1):
        result.append([result[i-1][0] + gap])
        for j in range(m_len + 1):
            break
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
