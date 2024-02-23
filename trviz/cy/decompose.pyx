import numpy as np
cimport numpy as np
cimport cython

# https://cython.readthedocs.io/en/latest/src/tutorial/numpy.html
np.import_array()

DTYPE = np.int32
ctypedef np.int32_t DTYPE_t
NEGATIVE_INF = -2147483648

# Define a Cython helper function to check if DP parameters are valid
cpdef void check_if_dp_parameters_are_valid(dict kwargs):
    cdef set valid_keys = {"match_score", "mismatch_score",
                           "insertion_score", "deletion_score", "min_score_threshold",
                           "verbose"}
    for k, v in kwargs.items():
        if k not in valid_keys:
            raise KeyError(f"Invalid keyword: {k}")
        if "score" in k:
            if type(v) not in [int, float]:
                raise ValueError(f"Invalid value {k}: {v}")
        if "verbose" == k:
            if type(v) != bool:
                raise ValueError(f"Invalid value type. {k}: {v}")

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef list decompose_cy(
        str sequence,
        list motifs,
        dict kwargs):

    check_if_dp_parameters_are_valid(kwargs)

    # Parameter setting
    cdef int match_score = kwargs.get("match_score", 1)
    cdef int mismatch_score = kwargs.get("mismatch_score", -1)
    cdef int insertion_score = kwargs.get("insertion_score", -1)
    cdef int deletion_score = kwargs.get("deletion_score", -1)
    cdef DTYPE_t min_score_threshold = kwargs.get("min_score_threshold", NEGATIVE_INF)
    cdef bint verbose = kwargs.get("verbose", False)

    # Define s[i, m, j]
    # Use C-style arrays instead of NumPy arrays
    # cdef int[:, :, :] s
    # cdef object[:, :, :] backtrack
    cdef np.ndarray[DTYPE_t, ndim=3] s
    cdef np.ndarray[object, ndim=3] backtrack
    cdef int max_motif_length = len(max(motifs, key=len))
    if verbose:
        print("Motifs used for decomposition: {}".format(','.join(motifs)))
        print("Max motif length", max_motif_length)

    # Allocate memory for s and backtrack arrays
    s = np.zeros((len(sequence) + 1, len(motifs), max_motif_length + 1), dtype=DTYPE)
    backtrack = np.zeros((len(sequence) + 1, len(motifs), max_motif_length + 1), dtype=object)

    # Boundary cases, when i = 0 or j = 0
    cdef int m, i, j
    for m, motif in enumerate(motifs):
        for i in range(len(sequence) + 1):
            for j in range(len(motif) + 1):
                if i == 0 and j == 0:
                    s[0, m, 0] = 0
                    backtrack[0, m, j] = (0, m, 0)
                elif i == 0 and j != 0:
                    s[0, m, j] = s[0, m, j - 1] + insertion_score
                    backtrack[0, m, j] = (0, m, j - 1)
                elif i != 0 and j == 0:
                    s[i, m, 0] = s[i - 1, m, 0] + insertion_score
                    backtrack[i, m, 0] = (i - 1, m, 0)

    # Normal cases, if i != 0
    # Note that sequence and motif are 1-based, but the DP table is 0-based.
    cdef DTYPE_t diagonal, from_diagonal, from_m_left, from_m_up, max_from_motif_end, max_motif_val
    cdef int argmax_index, max_m_index, max_j_of_max_m
    cdef int mi
    cdef str ms

    for i in range(1, len(sequence) + 1):
        for m, motif in enumerate(motifs):
            for j in range(1, len(motif) + 1):
                if j == 1:
                    if i == 1:
                        from_diagonal = s[i - 1, m, j - 1] + match_score if sequence[i - 1] == motif[j - 1] else \
                            s[i - 1, m, j - 1] + mismatch_score
                        from_m_left = s[i - 1, m, j] + insertion_score
                        from_m_up = s[i, m, j - 1] + deletion_score

                        # s[1][m][1]
                        s[i, m, j] = max(from_diagonal, from_m_left, from_m_up)
                        argmax_index = np.argmax([from_diagonal, from_m_left, from_m_up])

                        if argmax_index == 0:  # max from end
                            backtrack[i][m][j] = (0, m, 0)
                        elif argmax_index == 1:  # max from left in the same m
                            backtrack[i][m][j] = (i - 1, m, j)
                        else:  # max from up in the same m
                            backtrack[i][m][j] = (i, m, 0)  # j == 0
                    else:
                        max_motif_val = NEGATIVE_INF
                        max_m_index = -1
                        max_j_of_max_m = -1
                        for mi, ms in enumerate(motifs):
                            m_end = s[i - 1, mi, len(ms)]
                            if m_end > max_motif_val:
                                max_motif_val = m_end
                                max_m_index = mi
                                max_j_of_max_m = len(ms)

                        max_from_motif_end = max_motif_val + match_score if sequence[i - 1] == motif[
                            0] else max_motif_val + mismatch_score
                        from_m_left = s[i - 1, m, 1] + insertion_score
                        from_m_up = s[i, m, 0] + deletion_score

                        s[i, m, j] = max(max_from_motif_end, from_m_left, from_m_up)
                        '''
                            Tie-breaking for max_from_motif_end vs from_m_left
                            np.argmax select the first item if they have the same values
                            Select motif_end when motif_end and from_m_left has the same value
                            
                            In the following example, 2) is preferred.
                            Motif: ACGT, Sequence: ACGT"A"ACGT (Additional T in between the motif)
                            Decomposition 1) ACGT- AACGT
                            Decomposition 2) ACGTA -ACGT
                        '''
                        argmax_index = np.argmax([max_from_motif_end, from_m_left, from_m_up])

                        if argmax_index == 0:  # max from motif end
                            backtrack[i, m, j] = (i - 1, max_m_index, max_j_of_max_m)
                        elif argmax_index == 1:  # max from left in the same m
                            backtrack[i, m, j] = (i - 1, m, 1)
                        else:  # max from up in the same m
                            backtrack[i, m, j] = (i, m, 0)
                else:
                    diagonal = s[i - 1, m, j - 1] + match_score if sequence[i - 1] == motif[j - 1] else \
                    s[i - 1, m, j - 1] + mismatch_score
                    from_left = s[i - 1, m, j] + insertion_score
                    from_up = s[i, m, j - 1] + deletion_score

                    s[i, m, j] = max(diagonal, from_left, from_up)
                    argmax_index = np.argmax([diagonal, from_left, from_up])
                    if argmax_index == 0:
                        backtrack[i, m, j] = (i - 1, m, j - 1)
                    elif argmax_index == 1:
                        backtrack[i, m, j] = (i - 1, m, j)
                    else:
                        backtrack[i, m, j] = (i, m, j - 1)

    if verbose:
        print("DP table")
        for i in range(len(sequence) + 1):
            for m in range(len(motifs)):
                print(f"i{i}, m{m}, {s[i]}")

        print("Backtrack table")
        for i in range(len(sequence) + 1):
            for m in range(len(motifs) + 1):
                print(f"i{i}, m{m}, {backtrack[i]}")

    # Backtracking - getting decomposed motifs
    cdef object backtrack_start = None
    cdef DTYPE_t backtrack_max = min_score_threshold
    for m, motif in enumerate(motifs):
        if backtrack_max < s[len(sequence), m, len(motif)]:
            backtrack_max = s[len(sequence), m, len(motif)]
            backtrack_start = (len(sequence), m, len(motif))

    if backtrack_start is None:
        raise ValueError("No good match greater than score threshold of {}".format(min_score_threshold))

    if verbose:
        print("Best score: ", backtrack_max)

    cdef object backtrack_pointer = backtrack_start
    cdef int prev_i = -1
    cdef int prev_j = -1
    cdef int motif_end = len(sequence)
    cdef str decomposed_motif = ""
    cdef list decomposed_motifs = []

    while True:
        if verbose:
            print("Backtrack pointer", backtrack_pointer)
        i, m, j = backtrack_pointer

        if prev_j == 1 and j != 1:  # decompose
            # This condition indicates there was a jump. So, we need to decompose the motif
            # sequence[prev_i: motif_end].
            # update motif_end
            if verbose:
                print("Decomposed motif: ", decomposed_motif[::-1])
            if prev_i == 0:
                decomposed_motifs.append(sequence[:motif_end])
                # print("Added {}".format(sequence[:motif_end]))
            else:
                decomposed_motifs.append(sequence[prev_i-1:motif_end])
                # print("Added {}".format(sequence[prev_i-1:motif_end]))
            motif_end = prev_i-1 if prev_i > 0 else motif_end
            # decomposed_motifs.append(decomposed_motif[::-1])
            # decomposed_motif = ""

        # if prev_i != i and i != 0:  # emitted a symbol and not the start
        #     decomposed_motif += sequence[i - 1]

        if i == 0 and j == 0:
            # print("Backtrack finished")
            # print("motif_end", motif_end)
            # if motif_end > 1:
            #     print("Added {}".format(sequence[:motif_end]))
            #     decomposed_motifs.append(sequence[:motif_end])
            break

        prev_i = i
        prev_j = j
        backtrack_pointer = backtrack[i, m, j]

    if verbose:
        print("Input     : {}".format(''.join(decomposed_motifs[::-1])))
        print("Decomposed: {}".format(' '.join(decomposed_motifs[::-1])))

    return decomposed_motifs[::-1]
