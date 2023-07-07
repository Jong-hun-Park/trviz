import string
from collections import Counter

import numpy as np
import itertools
from Bio import SeqIO

LOWERCASE_LETTERS = string.ascii_lowercase
UPPERCASE_LETTERS = string.ascii_uppercase
DIGITS = string.digits

skipping_characters = ['(', '=', '<', '>', '?', '-']
PRIVATE_MOTIF_LABEL = '?'
INDEX_TO_CHR = list(LOWERCASE_LETTERS) + list(UPPERCASE_LETTERS) + list(DIGITS)
INDEX_TO_CHR.extend([chr(x) for x in range(33, 127) if chr(x) not in skipping_characters and chr(x) not in INDEX_TO_CHR])

DNA_CHARACTERS = {'A', 'C', 'G', 'T'}


def get_sample_and_sequence_from_fasta(fasta_file):
    """ Read fasta file and output headers and sequences """
    headers = []
    sequences = []
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            headers.append(record.id)
            sequences.append(str(record.seq.upper()))

    return headers, sequences

def get_motif_counter(decomposed_vntrs):
    """ Return a counter for each motif """
    motif_counter = Counter()
    for decomposed_vntr in decomposed_vntrs:
        motif_counter.update(Counter(decomposed_vntr))

    return motif_counter


def is_emitting_state(state_name):
    """ Check if the given state is emitting state, that is insertion or matching state """
    if state_name.startswith('M') or state_name.startswith('I') or state_name.startswith('start_random_matches') \
            or state_name.startswith('end_random_matches'):
        return True
    return False


def get_repeating_pattern_lengths(visited_states):
    lengths = []
    prev_start = None
    for i in range(len(visited_states)):
        if visited_states[i].startswith('unit_end') and prev_start is not None:
            current_len = 0
            for j in range(prev_start, i):
                if is_emitting_state(visited_states[j]):
                    current_len += 1
            lengths.append(current_len)
        if visited_states[i].startswith('unit_start'):
            prev_start = i
    return lengths


def get_motifs_from_visited_states_and_region(visited_states, region):
    lengths = get_repeating_pattern_lengths(visited_states)
    repeat_segments = []
    added = 0
    for l in lengths:
        repeat_segments.append(region[added:added + l])
        added += l
    return repeat_segments


def is_valid_sequence(sequence):
    """ Check if the given sequence is DNA sequence """
    for s in sequence:
        if s not in DNA_CHARACTERS:
            return False
    return True


def sort_by_manually(aligned_vntrs, sample_ids, sample_order_file):
    """ Sort the aligned and encoded tandem repeats based on the given order """
    with open(sample_order_file) as f:
        sample_order = [line.strip() for line in f.readlines()]

    sorted_sample_ids = []
    sorted_aligned_vntrs = []
    for sample_id in sample_order:
        if sample_id in sample_ids:
            sorted_sample_ids.append(sample_id)
            sorted_aligned_vntrs.append(aligned_vntrs[sample_ids.index(sample_id)])

    return sorted_sample_ids, sorted_aligned_vntrs


def sort(aligned_vntrs, sample_ids, symbol_to_motif, sample_order_file, method='motif_count'):
    """ Sort the aligned and encoded tandem repeats """
    if method == 'name':
        return zip(*sorted(list(zip(sample_ids, aligned_vntrs)), key=lambda x: x[0]))
    elif method == 'motif_count':
        return zip(*sorted(list(zip(sample_ids, aligned_vntrs)), key=lambda x: len(x[1].replace('-', ''))))
    elif method == 'simulated_annealing':
        return sort_by_simulated_annealing_optimized(aligned_vntrs, sample_ids, symbol_to_motif)
    elif method == 'manually':
        return sort_by_manually(aligned_vntrs, sample_ids, sample_order_file)
    else:
        raise ValueError("Please check the rearrangement method. {}".format(method))


def get_levenshtein_distance(s1, s2):
    """
    This function takes two strings and returns the Levenshtein distance between them.
    The Levenshtein distance is the minimum number of single-character edits (insertions, deletions or substitutions)
    required to change one string into the other.
    For example, the Levenshtein distance between "kitten" and "sitting" is 3, since the following three edits change
    one into the other, and there is no way to do it with fewer than three edits:
    kitten → sitten (substitution of "s" for "k")
    sitten → sittin (substitution of "i" for "e")
    """
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]


def _calculate_cost(seq1, seq2, alphabet_to_motif):
    if len(seq1) != len(seq2):
        raise Exception("The length of two sequences should be identical.")

    cost = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            # convert the motif to actual sequence
            if seq1[i] != '-' and seq2[i] != '-':
                s1 = alphabet_to_motif[seq1[i].lower()]
                s2 = alphabet_to_motif[seq2[i].lower()]
                cost += get_levenshtein_distance(s1, s2)
            else:
                if seq1[i] == '-':
                    cost += len(alphabet_to_motif[seq2[i].lower()])
                else:
                    cost += len(alphabet_to_motif[seq1[i].lower()])
    return cost


def calculate_cost_with_dist_matrix(aligned_encoded_vntr1, aligned_encoded_vntr2, dist_matrix, allow_copy_change=False):
    if len(aligned_encoded_vntr1) != len(aligned_encoded_vntr2):
        raise Exception("The length of two sequences should be identical.")
    cost = 0
    for i in range(len(aligned_encoded_vntr1)):
        symbol1 = aligned_encoded_vntr1[i]
        symbol2 = aligned_encoded_vntr2[i]
        if symbol1 != symbol2:
            # convert the motif to actual sequence
            if symbol1 != '-' and symbol2 != '-':
                cost += dist_matrix[symbol1][symbol2]
            else:
                if symbol1 == '-':
                    if allow_copy_change:
                        cost += 1  # allow motif count change as "one edit"
                    else:
                        cost += dist_matrix[symbol2][symbol2]  # length of encoded character 2
                else:
                    if allow_copy_change:
                        cost += 1  # allow motif count change as "one edit"
                    else:
                        cost += dist_matrix[symbol1][symbol1]  # length of encoded character 1

    return cost


def calculate_cost(alinged_vntrs, alphabet_to_motif):
    total_cost = 0
    for i in range(len(alinged_vntrs) - 1):
        total_cost += _calculate_cost(alinged_vntrs[i], alinged_vntrs[i + 1], alphabet_to_motif)

    return total_cost


def get_distance_matrix(symbol_to_motif, score=False):
    """
    Stores the edit distance between a motif and another motif.
    if two motifs are the same (e.g. dist_matrix[motif_x][motif_x]) it stores the length of the motif.

    :param symbol_to_motif: a dictionary mapping symbols to motifs
    :param score: if True, it outputs score matrix (1 - distance/max_dist)
    """
    dist_matrix = dict()
    max_score = 5

    for symbol1 in symbol_to_motif:
        dist_matrix[symbol1] = dict()
        for symbol2 in symbol_to_motif:
            motif_seq1 = symbol_to_motif[symbol1]
            motif_seq2 = symbol_to_motif[symbol2]
            if symbol1 == symbol2:
                dist_matrix[symbol1][symbol2] = len(motif_seq1)
            else:
                edit_dist = get_levenshtein_distance(motif_seq1, motif_seq2)
                if score:
                    dist_matrix[symbol1][symbol2] = max_score * (1.0 - edit_dist / max(len(motif_seq1), len(motif_seq2)))
                else:
                    dist_matrix[symbol1][symbol2] = edit_dist

    return dist_matrix


def get_score_matrix(symbol_to_motif,
                     match_score=2,
                     mismatch_score_for_edit_dist_of_1=-1,
                     mismatch_score_for_edit_dist_greater_than_1=-2,
                     gap_open_penalty=1.5,
                     gap_extension_penalty=0.6,
                     ):

    score_matrix = dict()
    score_matrix['gap_open'] = gap_open_penalty
    score_matrix['gap_extension'] = gap_extension_penalty

    for symbol1 in symbol_to_motif:
        score_matrix[symbol1] = dict()
        for symbol2 in symbol_to_motif:
            motif_seq1 = symbol_to_motif[symbol1]
            motif_seq2 = symbol_to_motif[symbol2]
            if symbol1 == symbol2:
                score_matrix[symbol1][symbol2] = match_score
            else:
                edit_dist = get_levenshtein_distance(motif_seq1, motif_seq2)

                edit_dist_cutoff = 1
                if abs(len(motif_seq1) - len(motif_seq2)) <= 1:
                    edit_dist_cutoff += len(max(motif_seq2, motif_seq1, key=len)) // 30
                if edit_dist <= edit_dist_cutoff:
                    score_matrix[symbol1][symbol2] = mismatch_score_for_edit_dist_of_1
                else:
                    score_matrix[symbol1][symbol2] = mismatch_score_for_edit_dist_greater_than_1

    return score_matrix

def calculate_total_cost(alinged_vntrs, dist_matrix):
    total_cost = 0
    for i in range(len(alinged_vntrs) - 1):
        total_cost += calculate_cost_with_dist_matrix(alinged_vntrs[i], alinged_vntrs[i + 1], dist_matrix)

    return total_cost


def sort_by_simulated_annealing_optimized(seq_list, sample_ids, symbol_to_motif):
    dist_matrix = get_distance_matrix(symbol_to_motif)

    initial_cost = calculate_total_cost(seq_list, dist_matrix)
    initial_seq_list = seq_list.copy()
    initial_sample_ids = sample_ids.copy()
    initial_cost_test = calculate_cost(seq_list, symbol_to_motif)
    # assert initial_cost == initial_cost_test, "Should be the same {} {}".format(initial_cost, initial_cost_test)

    T = 1_000
    DECAY = 0.9

    iteration = 0

    all_index_pairs = list(itertools.combinations(range(len(seq_list)), 2))
    while True:
        iteration += 1
        print("T:", T)
        if T <= 1e-2:
            break
        print("iteration", iteration)
        not_changed_count = 0

        # from random import choice
        # index_1, index_2 = choice(all_index_pairs)

        # from random import shuffle
        # shuffle(all_index_pairs)

        for index_1, index_2 in all_index_pairs:
            # only compare the cost before and after changing the order
            current_cost = 0
            after_cost = 0

            # Flanking cost for the index_1 sequence
            # Right side
            current_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_1 + 1], dist_matrix)
            if index_1 + 1 == index_2:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_1], dist_matrix)
            else:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_1 + 1], dist_matrix)
            if index_1 != 0:  # has left side
                # Left side
                current_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_1 - 1], dist_matrix)
                # index_1 < index_2, so index_1 -1 != index_2
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_1 - 1], dist_matrix)

            # Flanking cost for the index_2 sequence
            # Left side
            current_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_2 - 1], dist_matrix)
            if index_2 - 1 == index_1:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_2], dist_matrix)
            else:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_2 - 1], dist_matrix)
            if index_2 != len(seq_list) - 1:
                # Right side
                current_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_2 + 1], dist_matrix)
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_2 + 1], dist_matrix)

            if after_cost == current_cost:
                continue
            elif after_cost < current_cost:
                print("Swap occurred at {} and {}".format(index_1, index_2), "after cost", after_cost, "cur cost", current_cost)
                seq_list[index_1], seq_list[index_2] = seq_list[index_2], seq_list[index_1]
                sample_ids[index_1], sample_ids[index_2] = sample_ids[index_2], sample_ids[index_1]
            else:
                prob = np.exp(-(after_cost - current_cost) * 10 / T)
                if prob > np.random.uniform(low=0.0, high=1.0):
                    # swap
                    print("Swap occurred {}".format(prob), "after cost", after_cost, "cur cost", current_cost)
                    seq_list[index_1], seq_list[index_2] = seq_list[index_2], seq_list[index_1]
                    sample_ids[index_1], sample_ids[index_2] = sample_ids[index_2], sample_ids[index_1]
                else:
                    not_changed_count += 1

        T *= DECAY
        # if not_changed_count == len(seq_list):
        #     break

    print("The initial cost", initial_cost)
    after_cost = calculate_total_cost(seq_list, dist_matrix)
    print("Cost after sorting", after_cost)

    if initial_cost < after_cost:
        # print("The result has higher cost, so it doesn't change the order")
        return initial_sample_ids, initial_seq_list
    return sample_ids, seq_list


def add_padding(encoded_trs):
    """
    This function takes a list of encoded traces as input and returns a list of padded traces.
    The padding is done by adding '-' to the end of each trace.
    The number of '-' added to each trace is equal to the difference between the length of the longest trace and
    the length of the trace.
    """
    max_motif_count = len(max(encoded_trs, key=len))
    padded_trs = []
    for encoded_tr in encoded_trs:
        padding_count = max_motif_count - len(encoded_tr)
        padded_trs.append(encoded_tr + '-' * padding_count)

    return padded_trs


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', print_end = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=print_end)
    # Print New Line on Complete
    if iteration == total:
        print()
