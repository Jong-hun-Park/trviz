from collections import Counter

from trviz.utils import is_valid_sequence
from trviz.utils import get_motifs_from_visited_states_and_region

from trviz.utils import INDEX_TO_CHR

import settings
if settings.decomposition_method == "hmm":
    from pomegranate import DiscreteDistribution, State
    from pomegranate import HiddenMarkovModel as Model

import numpy as np


class TandemRepeatDecomposer:

    def __init__(self, mode="DP"):
        if mode == "DP" or mode == "HMM":
            self.mode = mode
        else:
            raise ValueError(f"{mode} is invalid mode for tandem repeat decomposer.")

    # def decompose(self):
    #     if self.mode == "DP":
    #         self.decompose_dp()
    #     else:
    #         self.decompose_hmm()

    @staticmethod
    def _build_repeat_finder_hmm(motif, copies=1, has_flanking_sequence=False):
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

    # Problem: Tandem repeat motif labeling
    # Input:
    # 1. tandem repeat sequences with the begin and end of the mark. (two alleles for humans) with sample ID?
    # 2. (Optional) genomic context (flanking regions) of the TR sequences.
    # 3. (Optional) motif or motifs
    # Output:
    # 1. Set of motifs with numbers (labels)
    # 2. Decomposition of the string into the numbers (labels)

    # Input 1 is solved.
    # When 2 is given, we can use different types of HMM
    # When 3 is not given, we can find a consensus motif using TRF.
    @staticmethod
    def decompose_dp(sequence,
                     motifs,
                     match_score=5,
                     mismatch_score=-2,
                     min_score_threshold=float("-inf"),
                     verbose=False):
        """
        Decompose sequence into motifs using dynamic programming
        :param sequence: a string of VNTR sequence
        :param motifs: a list of motif strings composed of nucleotide strings

        :return: A list of decomposed motifs
        """

        if not is_valid_sequence(sequence):
            raise ValueError("Sequence has invalid characters")

        if isinstance(motifs, str):
            motifs = [motifs]  # only one string is given

        for motif in motifs:
            if not is_valid_sequence(motif):
                raise ValueError("Consensus motif has invalid characters")

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
            -  
        """

        # Score
        insertion_score = mismatch_score
        deletion_score = mismatch_score

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

        s = np.zeros(len(sequence) + 1, dtype=object)
        backtrack = np.zeros(len(sequence) + 1, dtype=object)
        for i in range(len(sequence) + 1):
            s[i] = np.zeros(len(motifs), dtype=object)
            backtrack[i] = np.zeros(len(motifs), dtype=object)
            for m in range(len(motifs)):
                s[i][m] = np.zeros(max_motif_length + 1)
                backtrack[i][m] = np.zeros(max_motif_length + 1, dtype=object)

        # Boundary cases, when i = 0 or j = 0
        for m, motif in enumerate(motifs):
            for i in range(len(sequence) + 1):
                for j in range(len(motif) + 1):
                    if i == 0 and j == 0:
                        s[0][m][0] = 0
                        backtrack[0][m][j] = (0, m, 0)
                    elif i == 0 and j != 0:
                        s[0][m][j] = s[0][m][j - 1] + insertion_score
                        backtrack[0][m][j] = (0, m, j - 1)
                    elif i != 0 and j == 0:
                        s[i][m][0] = s[i - 1][m][0] + insertion_score
                        backtrack[i][m][0] = (i - 1, m, 0)

        # Normal cases, if i != 0
        for i in range(1, len(sequence) + 1):
            for m, motif in enumerate(motifs):
                for j in range(1, len(motif) + 1):
                    if j == 1:
                        if i == 1:
                            from_diagonal = s[i - 1][m][j - 1] + match_score if sequence[i - 1] == motif[j - 1] else \
                            s[i - 1][m][j - 1] + mismatch_score
                            from_m_left = s[i - 1][m][1] + match_score if sequence[i - 1] == motif[j - 1] else \
                            s[i - 1][m][1] + mismatch_score
                            from_m_up = s[i][m][j - 1] + match_score if sequence[i - 1] == motif[0] else s[i][m][
                                                                                                             j - 1] + mismatch_score

                            s[i][m][j] = max(from_diagonal, from_m_left, from_m_up)
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
                                m_end = s[i - 1][mi][len(ms)]
                                if m_end > max_motif_val:
                                    max_motif_val = m_end
                                    max_m_index = mi
                                    max_j_of_max_m = len(ms)

                            max_from_end = max_motif_val + match_score if sequence[i - 1] == motif[
                                0] else max_motif_val + mismatch_score
                            from_m_left = s[i - 1][m][1] + match_score if sequence[i - 1] == motif[0] else s[i - 1][m][1] + mismatch_score
                            from_m_up = s[i][m][0] + match_score if sequence[i - 1] == motif[0] else s[i][m][0] + mismatch_score

                            s[i][m][j] = max(max_from_end, from_m_left, from_m_up)
                            argmax_index = np.argmax([max_from_end, from_m_left, from_m_up])
                            if argmax_index == 0:  # max from end
                                backtrack[i][m][j] = (i - 1, max_m_index, max_j_of_max_m)
                            elif argmax_index == 1:  # max from left in the same m
                                backtrack[i][m][j] = (i - 1, m, 1)
                            else:  # max from up in the same m
                                backtrack[i][m][j] = (i, m, 0)

                    else:
                        # print(f'motif {motif}, i {i}, m {m}, j{j}')
                        diagonal = s[i - 1][m][j - 1] + match_score if sequence[i - 1] == motif[j - 1] else s[i - 1][m][j - 1] + mismatch_score
                        insertion = s[i - i][m][j] + insertion_score
                        deletion = s[i][m][j - 1] + deletion_score

                        s[i][m][j] = max(diagonal, insertion, deletion)
                        path = np.argmax([diagonal, insertion, deletion])
                        if path == 0:
                            backtrack[i][m][j] = (i - 1, m, j - 1)
                        elif path == 1:
                            backtrack[i][m][j] = (i - 1, m, j)
                        else:
                            backtrack[i][m][j] = (i, m, j - 1)

        if verbose:
            print("DP table")
            for i in range(len(sequence)):
                for m in range(len(motifs)):
                    print(f"i{i}, m{m}, {s[i]}")

            print("Backtrack table")
            # Print backtrack list
            for i in range(len(sequence)):
                for m in range(len(motifs)):
                    print(f"i{i}, m{m}, {backtrack[i]}")

        # Backtracking - getting decomposed motifs
        backtrack_start = None
        backtrack_max = min_score_threshold
        for m, motif in enumerate(motifs):
            if backtrack_max < s[len(sequence)][m][len(motif)]:
                backtrack_max = s[len(sequence)][m][len(motif)]
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

            if prev_i != i and i != 0:
                decomposed_motif += sequence[i - 1]

            backtrack_pointer = backtrack[i][m][j]

            if i == 0 and j == 0:
                break

            prev_i = i
            prev_j = j

        if verbose:
            print("Input     : {}".format(''.join(decomposed_motifs[::-1])))
            print("Decomposed: {}".format(' '.join(decomposed_motifs[::-1])))

        return decomposed_motifs[::-1]

    def decompose_hmm(self, sequence, consensus_motif=None, repeat_count=None, has_flanking_sequence=False, verbose=False):
        """
        Decompose sequence into motifs using a HMM
        :param sequence: a string of VNTR sequence
        :param consensus_motif: A motif sequence
        :param repeat_count: If not specified, use estimated repeat count based on the sequence and motif length
        :param has_flanking_sequence: If true, use additional nodes in the HMM to identify flanking regions
        :return: A list of decomposed motifs
        """

        if not is_valid_sequence(sequence):
            raise ValueError("Sequence has invalid characters")

        if consensus_motif is not None and not is_valid_sequence(consensus_motif):
            raise ValueError("Consensus motif has invalid characters")

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

    @staticmethod
    def get_motif_counter(decomposed_vntrs):
        motif_counter = Counter()
        for decomposed_vntr in decomposed_vntrs:
            motif_counter.update(Counter(decomposed_vntr))

        return motif_counter

    @staticmethod
    def _divide_motifs_into_normal_and_private(decomposed_vntrs, private_motif_threshold):
        """
        Givne a list of decomposed VNTRs, divide motifs into two groups: normal and private.
        If a motif occurred less than the private_motif_threshold, it is regarded as private motif.
        Otherwise, normal motifs.

        :param decomposed_vntrs
        :param private_motif_threshold
        :return: normal motifs, private motifs
        """

        motif_counter = TandemRepeatDecomposer.get_motif_counter(decomposed_vntrs)

        normal_motifs = []
        private_motifs = []
        for motif, count in motif_counter.items():
            if count > private_motif_threshold:
                normal_motifs.append(motif)
            else:
                private_motifs.append(motif)

        return normal_motifs, private_motifs, motif_counter

    @staticmethod
    def find_minimum_private_motif_threshold(decomposed_vntrs):
        motif_counter = TandemRepeatDecomposer.get_motif_counter(decomposed_vntrs)

        min_private_motif_threshold = 0
        for index, (motif, count) in enumerate(motif_counter.most_common()):
            print(index, motif, count)
            if index + 1 > len(INDEX_TO_CHR) - 1:  # starting with 33, last 126, skipping 4 symbols, 1 for private
                min_private_motif_threshold = count
                break

        return min_private_motif_threshold

    @staticmethod
    def label_motifs(decomposed_vntrs, private_motif_threshold=0, auto=None):
        """

        :param decomposed_vntrs:
        :param private_motif_threshold:
        :param auto: if True, find the minimum threshold to encode everything using 126-33+1 ASCII characters.
        :return: labeled_vntrs, motif to alphabet (dictionary of the mapping)
        """

        def _get_motif_labels(decomposed_vntrs):
            """
            :param decomposed_vntrs
            :return: unique motifs
            """
            unique_motifs = set()
            for vntr in decomposed_vntrs:
                for motif in vntr:
                    unique_motifs.add(motif)

            return unique_motifs

        def _index_to_char(index):
            """
            --anysymbol
            To use unusual characters (e.g., U as selenocysteine in protein sequence; i as inosine in nucleotide sequence),
            use the --anysymbol option:
            % mafft --anysymbol input > output
            It accepts any printable characters (U, O, #, $, %, etc.; 0x21-0x7e in the ASCII code),
            execpt for > (0x3e) and ( (0x28).
            # '= (60) < (61) > (62)' can not be used
            Unusual characters are scored as unknown (not considered in the calculation), unlike in the --text mode.
            """
            if index < 0 or index > len(INDEX_TO_CHR) - 1:
                raise ValueError(f"Index should range between 0 to {len(INDEX_TO_CHR) - 1}. Given : {index}")

            return INDEX_TO_CHR[index]

        if auto:
            private_motif_threshold = TandemRepeatDecomposer.find_minimum_private_motif_threshold(decomposed_vntrs)
            print("private motif threshold: ", private_motif_threshold)

        motif_to_symbol = {}
        symbol_to_motif = {}
        # For private motifs, we use single letter to encode them.
        if private_motif_threshold > 0:
            normal_motifs, private_motifs, motif_counter = TandemRepeatDecomposer._divide_motifs_into_normal_and_private(decomposed_vntrs,
                                                                                                 private_motif_threshold)
            if len(normal_motifs) + 1 > len(INDEX_TO_CHR):
                print("Motif counter:", motif_counter)
                raise ValueError("Too many unique motifs. Can not encode properly: {} unique motifs".format(
                    len(normal_motifs) + len(private_motifs)))

            # debug
            print(motif_counter)
            print(normal_motifs)
            print(private_motifs)

            # Assign a code to all private motifs
            motif_to_symbol.update({motif: _index_to_char(0) for motif in sorted(private_motifs)})
            symbol_to_motif.update({_index_to_char(0): motif for motif in sorted(private_motifs)})

            # For normal motifs
            motif_to_symbol.update(
                {motif: _index_to_char(index + 1) for index, motif in enumerate(sorted(normal_motifs))})
            symbol_to_motif.update(
                {_index_to_char(index + 1): motif for index, motif in enumerate(sorted(normal_motifs))})
        else:  # Use all distinct motif
            unique_motifs = _get_motif_labels(decomposed_vntrs)
            if len(unique_motifs) > len(INDEX_TO_CHR):  # single symbol ascii
                raise ValueError(
                    "Too many unique motifs. Can not encode properly: {} unique motifs".format(len(unique_motifs)))

            motif_to_symbol.update(
                {motif: _index_to_char(index) for index, motif in enumerate(sorted(list(unique_motifs)))})
            symbol_to_motif.update(
                {_index_to_char(index): motif for index, motif in enumerate(sorted(list(unique_motifs)))})

        # Label VNTRs using the encoding
        labeled_vntrs = []
        for vntr in decomposed_vntrs:
            labeled_vntr = ""
            for motif in vntr:
                labeled_vntr += str(motif_to_symbol[motif])
            labeled_vntrs.append(labeled_vntr)

        return labeled_vntrs, motif_to_symbol, symbol_to_motif


if __name__ == "__main__":
    pass
    # sequence = "ACTGGGACTGACTGT"
    # motifs = ["ACTG", "ACTGGG", "ACTGT"]
    # decomposed_motifs = decompose_dp(sequence, motifs, verbose=True)
    # print(f'sequence {sequence}')
    # print(f'motifs {motifs}')
    # print(f'decomposed_motifs {decomposed_motifs}')
    #
    # sequence = "AAATAAAATAAAATAAAATA"
    # #          "AAATA AAATTA AAATA AAAAATA"
    # # motifs = ["AAATA"]
    # motifs = "AAATA"
    # decomposed_motifs = decompose_dp(sequence, motifs, verbose=True)
    # print(f'sequence {sequence}')
    # print(f'motifs {motifs}')
    # print(f'decomposed_motifs {decomposed_motifs}')
    # exit(1)
    #
    # decomposed_vntrs = []
    #
    # sequence = "ACTGACTGACTG"
    # consensus_motif = "ACTG"
    # # decomposed_motifs = decompose_hmm(sequence, consensus_motif)
    # decomposed_motifs = decompose_dp(sequence, consensus_motif)
    # print(decomposed_motifs)
    # decomposed_vntrs.append(decomposed_motifs)
    # # Collect multiple decomposed motifs.
    #
    # sequence = "ACTGACTGACTGACCTGACTG"
    # consensus_motif = "ACTG"
    # # decomposed_motifs = decompose_hmm(sequence, consensus_motif)
    # decomposed_motifs = decompose_dp(sequence, consensus_motif)
    # print(decomposed_motifs)
    # decomposed_vntrs.append(decomposed_motifs)
    #
    # sequence = "ACTGACTGACTTGACCTGACTGACTGACTG"
    # consensus_motif = "ACTG"
    # # decomposed_motifs = decompose_hmm(sequence, consensus_motif)
    # decomposed_motifs = decompose_dp(sequence, consensus_motif)
    # print(decomposed_motifs)
    # decomposed_vntrs.append(decomposed_motifs)
    #
    # labeled_vntrs, motif_to_alphabet = label_motifs(decomposed_vntrs)
    # print(motif_to_alphabet)
    # print(labeled_vntrs)
    #
    # # Alignment test
    # from trviz.motif_aligner import align_motifs
    #
    # aligned_vntrs = align_motifs(labeled_vntrs)
    # print(aligned_vntrs)
    #
    # # Visualization test
    # from trviz.visualization import trplot
    #
    # trplot(aligned_vntrs)
