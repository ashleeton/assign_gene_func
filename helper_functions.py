from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")

GAP_PENALTY = -8

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    n = len(seq1)
    m = len(seq2)
    matrix = [[0] * (m + 1) for i in range(n + 1)]

    for i in range(1, n + 1):
        matrix[i][0] = matrix[i - 1][0] + scoring_function(seq1[i - 1], "-")

    for j in range(1, m + 1):
        matrix[0][j] = matrix[0][j - 1] + scoring_function(seq2[j - 1], "-")

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = matrix[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = matrix[i - 1][j] + scoring_function(seq1[i - 1], "-")
            down = matrix[i][j - 1] + scoring_function(seq2[j - 1], "-")
            matrix[i][j] = max(diagonal, up, down)

    aligned1 = ""
    aligned2 = ""
    i = n
    j = m

    while i > 0 or j > 0:
        if matrix[i][j] == matrix[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1]):
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i - 1][j] + scoring_function(seq1[i - 1], "-"):
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        elif matrix[i][j] == matrix[i][j - 1] + scoring_function(seq2[j - 1], "-"):
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1

    return aligned1, aligned2, float(H[n][m])


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    n = len(seq1)
    m = len(seq2)
    matrix = [[0] * (m + 1) for i in range(n + 1)]

    for i in range(1, n + 1):
        matrix[i][0] = 0

    for j in range(1, m + 1):
        matrix[0][j] = 0

    best_score = 0
    best_i = 0
    best_j = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = matrix[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = matrix[i - 1][j] + scoring_function(seq1[i - 1], "-")
            down = matrix[i][j - 1] + scoring_function(seq2[j - 1], "-")
            matrix[i][j] = max(0, diagonal, up, down)
            if matrix[i][j] > best_score:
                best_score = matrix[i][j]
                best_i = i
                best_j = j

    aligned1 = ""
    aligned2 = ""
    i = best_i
    j = best_j

    while matrix[i][j] != 0:
        if matrix[i][j] == matrix[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1]):
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i - 1][j] + scoring_function(seq1[i - 1], "-"):
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        elif matrix[i][j] == matrix[i][j - 1] + scoring_function(seq2[j - 1], "-"):
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1

    return aligned1, aligned2, float(best_score)


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

def scoring_function_blosum62(aa_i, aa_j):
    if aa_i == "-" or aa_j == "-":
        return GAP_PENALTY
    else:
        return blosum62[aa_i][aa_j]
