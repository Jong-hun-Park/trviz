from typing import List

from trviz.utils import is_valid_sequence
from trviz.utils import get_motifs_from_visited_states_and_region

import numpy as np


class Decomposer:

    def __init__(self, mode="DP"):
        if mode == "DP":
            self.mode = mode
        elif mode == "HMM":
            self.mode = mode
        else:
            raise ValueError(f"{mode} is invalid mode for tandem repeat decomposer.")

    def decompose(self, sequence, motifs, **kwargs):
        """
        Decompose sequence into motifs
        :param sequence: a string of tandem repeat sequence
        :param motifs: a set of motifs to be used for decomposition
        :param kwargs: keyword arguments for decomposer
        :return: decomposed sequence
        """
        if not isinstance(sequence, str):
            TypeError("Sequence must be a string")
        if isinstance(motifs, str):
            motifs = [motifs]  # only one string is given
        if not isinstance(motifs, list):
            TypeError("Motifs must be a list of strings")

        sequence = sequence.upper()
        motifs = [m.upper() for m in motifs]

        if not is_valid_sequence(sequence):
            raise ValueError(f"Sequence has invalid characters: {sequence}")

        for motif in motifs:
            if not is_valid_sequence(motif):
                raise ValueError(f"The motif has invalid characters: {motif}")

        if self.mode == "DP":
            return self._decompose_dp(sequence, motifs, **kwargs)
        else:
            return self._decompose_hmm(sequence, motifs, **kwargs)

    @staticmethod
    def _decompose_dp(
            sequence,
            motifs,
            **kwargs,
    ) -> List:
        """
        Decompose sequence into motifs using dynamic programming

        :param sequence: a string of VNTR sequence
        :param motifs: a list of motif strings composed of nucleotide strings
        :param **kwargs:
        valid keywords are the following:
        match_score: a score for match
        mismatch_score: a score for mismatch
        min_score_threshold: a minimum score of the alignment.

        :return: A list of decomposed motifs
        """

        """
        Let s[i,m,j] be the best parse of the sequence prefix s[1..i],
        where the i-th character is matched to the j-th character of motif m.
        Interesting case is when j=1 (1-st index of the motif)

        # Recurrence formulation

        # Boundary case
        For all m,
            if j = 0,
                s[0][m][0] = match_score if sequence[0] == motif[0] else mismatch_score
            if j != 0,
                s[0][m][j] = s[0][m][j-1] + insertion_score

        s[i,m,j] =
        1. if j = 1
                   max(max_m'( s[i-1, m', len(m')] ) + match )
                      (s[i-1, m, 1] + indel               )
        2. otherwise
                      ( s[i-1, m, j-1] + match    )
                   max( s[i-1, m, j]   + deletion )
                      ( s[i, m, j-1]   + insertion )


        # Solution
        The best parsing score: max_m( s[len(sequence), m, len(m)] )

        # Proof
        - i = 1, by definition
        - let i >= 2, j >= 2, assume that s[i][m][j] have been computed correctly.
        - We are claiming that the formulation for s[i, m, j] is the best answer
        """

        def _check_if_dp_parameters_are_valid(**kwargs):
            valid_keys = {"match_score", "mismatch_score",
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

        _check_if_dp_parameters_are_valid(**kwargs)

        # Parameter setting
        match_score = kwargs.get("match_score", 5)
        mismatch_score = kwargs.get("mismatch_score", -4)
        insertion_score = kwargs.get("insertion_score", -4)
        deletion_score = kwargs.get("deletion_score", -4)
        min_score_threshold = kwargs.get("min_score_threshold", float("-inf"))
        verbose = kwargs.get("verbose", False)

        # Define s[i, m, j]
        # numpy array doesn't allow to have different length array
        # Use a largest length of motifs to build a 3D square-like array
        # 0 <= i <= n: count: n+1
        # 0 <= m <= len(motif) - 1: count: len(motif)
        # 0 <= j <= len(m): count: len(m)+1
        max_motif_length = len(max(motifs, key=len))
        if verbose:
            print("Motifs used for decomposition: {}".format(','.join(motifs)))
            print("Max motif length", max_motif_length)

        s = np.zeros((len(sequence) + 1, len(motifs), max_motif_length + 1))
        backtrack = np.zeros((len(sequence) + 1, len(motifs), max_motif_length + 1), dtype=object)

        # Boundary cases, when i = 0 or j = 0
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
        for i in range(1, len(sequence) + 1):
            for m, motif in enumerate(motifs):
                for j in range(1, len(motif) + 1):
                    # print(f'motif: {motif}, i:{i}, m:{m}, j:{j}')
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
                            max_motif_val = float('-inf')
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
                                Motif: ACGT, Sequence: ACGT"T"ACGT (Additional T in between the motif)
                                Decomposition 1) ACGT- TACGT
                                Decomposition 2) ACGTT -ACGT
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
                        path = np.argmax([diagonal, from_left, from_up])
                        if path == 0:
                            backtrack[i, m, j] = (i - 1, m, j - 1)
                        elif path == 1:
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
        backtrack_start = None
        backtrack_max = min_score_threshold
        for m, motif in enumerate(motifs):
            if backtrack_max < s[len(sequence), m, len(motif)]:
                backtrack_max = s[len(sequence), m, len(motif)]
                backtrack_start = (len(sequence), m, len(motif))

        if backtrack_start is None:
            raise ValueError("No good match greater than score threshold of {}".format(min_score_threshold))

        if verbose:
            print("Best score: ", backtrack_max)

        backtrack_pointer = backtrack_start
        prev_i = -1
        prev_j = -1
        decomposed_motif = ""
        decomposed_motifs = []

        while True:
            if verbose:
                print("Backtrack pointer", backtrack_pointer)
            i, m, j = backtrack_pointer

            if prev_j == 1 and j != 1:  # decompose
                if verbose:
                    print("Decomposed motif: ", decomposed_motif[::-1])
                decomposed_motifs.append(decomposed_motif[::-1])
                decomposed_motif = ""
            if prev_i != i and i != 0:  # emitted a symbol and not the start
                decomposed_motif += sequence[i - 1]

            if i == 0 and j == 0:
                break
            backtrack_pointer = backtrack[i, m, j]
            prev_i = i
            prev_j = j

        if verbose:
            print("Input     : {}".format(''.join(decomposed_motifs[::-1])))
            print("Decomposed: {}".format(' '.join(decomposed_motifs[::-1])))

        return decomposed_motifs[::-1]

    @staticmethod
    def _build_repeat_finder_hmm(motif, copies=1, has_flanking_sequence=False):
        from pomegranate import DiscreteDistribution, State
        from pomegranate import HiddenMarkovModel as Model

        model = Model(name="RepeatFinderHMM")

        insert_distribution = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        if has_flanking_sequence:
            start_random_matches = State(insert_distribution, name='start_random_matches')
            end_random_matches = State(insert_distribution, name='end_random_matches')
            model.add_states([start_random_matches, end_random_matches])

        last_end = None
        for repeat in range(1, copies + 1):  # 1-based index
            insert_states = []
            match_states = []
            delete_states = []
            for i in range(len(motif) + 1):  # I0 to IN
                insert_states.append(State(insert_distribution, name='I%s_%s' % (i, repeat)))

            for i in range(len(motif)):  # M1 to MN
                distribution_map = dict({'A': 0.01, 'C': 0.01, 'G': 0.01, 'T': 0.01})
                distribution_map[motif[i]] = 0.97
                match_states.append(State(DiscreteDistribution(distribution_map), name='M%s_%s' % (str(i + 1), repeat)))

            for i in range(len(motif)):  # D1 to DN
                delete_states.append(State(None, name='D%s_%s' % (str(i + 1), repeat)))

            unit_start = State(None, name='unit_start_%s' % repeat)
            unit_end = State(None, name='unit_end_%s' % repeat)
            model.add_states(insert_states + match_states + delete_states + [unit_start, unit_end])
            last = len(delete_states) - 1

            if repeat == 1:  # The first repeat
                # If the sequence contains flanking sequences, we use random match states to distinguish them from TR seq
                if has_flanking_sequence:
                    model.add_transition(model.start, unit_start, 0.5)
                    model.add_transition(model.start, start_random_matches, 0.5)
                    model.add_transition(start_random_matches, unit_start, 0.5)
                    model.add_transition(start_random_matches, start_random_matches, 0.5)
                else:
                    model.add_transition(model.start, unit_start, 1)
            else:  # connect the previous unit end to unit start
                if has_flanking_sequence:
                    model.add_transition(last_end, unit_start, 0.5)
                else:
                    model.add_transition(last_end, unit_start, 1)

            if has_flanking_sequence:
                model.add_transition(unit_end, end_random_matches, 0.5)

            if repeat == copies:  # The last repeat
                if has_flanking_sequence:
                    model.add_transition(unit_end, model.end, 0.5)
                    model.add_transition(end_random_matches, end_random_matches, 0.5)
                    model.add_transition(end_random_matches, model.end, 0.5)
                else:
                    model.add_transition(unit_end, model.end, 1)

            model.add_transition(unit_start, match_states[0], 0.98)
            model.add_transition(unit_start, delete_states[0], 0.01)
            model.add_transition(unit_start, insert_states[0], 0.01)

            model.add_transition(insert_states[0], insert_states[0], 0.01)
            model.add_transition(insert_states[0], delete_states[0], 0.01)
            model.add_transition(insert_states[0], match_states[0], 0.98)

            model.add_transition(delete_states[last], unit_end, 0.99)
            model.add_transition(delete_states[last], insert_states[last + 1], 0.01)

            model.add_transition(match_states[last], unit_end, 0.99)
            model.add_transition(match_states[last], insert_states[last + 1], 0.01)

            model.add_transition(insert_states[last + 1], insert_states[last + 1], 0.01)
            model.add_transition(insert_states[last + 1], unit_end, 0.99)

            for i in range(0, len(motif)):
                model.add_transition(match_states[i], insert_states[i + 1], 0.01)
                model.add_transition(delete_states[i], insert_states[i + 1], 0.01)
                model.add_transition(insert_states[i + 1], insert_states[i + 1], 0.01)
                if i < len(motif) - 1:
                    model.add_transition(insert_states[i + 1], match_states[i + 1], 0.98)
                    model.add_transition(insert_states[i + 1], delete_states[i + 1], 0.01)

                    model.add_transition(match_states[i], match_states[i + 1], 0.98)
                    model.add_transition(match_states[i], delete_states[i + 1], 0.01)

                    model.add_transition(delete_states[i], delete_states[i + 1], 0.01)
                    model.add_transition(delete_states[i], match_states[i + 1], 0.98)

            last_end = unit_end

        model.bake(merge="None")  # Not merging silent states

        return model

    def _decompose_hmm(self, sequence, consensus_motif, **kwargs):
        """
        Decompose sequence into motifs using a HMM
        :param sequence: a string of VNTR sequence
        :param consensus_motif: A motif sequence
        :param **kwargs:

        1. repeat_count: If not specified, use estimated repeat count based on the sequence and motif length
        2. has_flanking_sequence: If true, use additional nodes in the HMM to identify flanking regions
        3. verbose

        :return: A list of decomposed motifs
        """

        def _check_if_hmm_parameters_are_valid(**kwargs):
            valid_keys = {"repeat_count", "has_flanking_sequence", "verbose"}
            for k, v in kwargs.items():
                if k not in valid_keys:
                    raise KeyError(f"Invalid keyword: {k}")
                if "repeat_count" in k:
                    if type(v) != int:
                        raise ValueError(f"Invalid value {k}: {v}")
                if "has_flanking_sequence" in k or "verbose" == k:
                    if type(v) != bool:
                        raise ValueError(f"Invalid value {k}: {v}")

        _check_if_hmm_parameters_are_valid(**kwargs)

        # Parameter setting
        repeat_count = kwargs.get("repeat_count", None)
        has_flanking_sequence = kwargs.get("has_flanking_sequence", False)
        verbose = kwargs.get("verbose", False)

        # if consensus_motif is None:
        #     #TODO Find a motif using TRF
        #     consensus_motif = find_motif_using_TRF(sequence)

        if repeat_count is None:
            repeat_count = round(len(sequence) / len(consensus_motif))

        if verbose:
            print("Estimated repeat count", repeat_count)
            print("Building HMM...")
            print("motif", consensus_motif)
            print("sequence", sequence)

        repeat_finder_hmm = self._build_repeat_finder_hmm(consensus_motif, repeat_count, has_flanking_sequence)
        if verbose:
            print("Building HMM done")
        logp, path = repeat_finder_hmm.viterbi(sequence)
        visited_states = [state.name for idx, state in path[1:-1]]
        deomposed_motifs = get_motifs_from_visited_states_and_region(visited_states, sequence)
        if verbose:
            print(logp)
            print(deomposed_motifs)

        return deomposed_motifs

