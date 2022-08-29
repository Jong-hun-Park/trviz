import string
from collections import Counter

import numpy as np
import itertools

LOWERCASE_LETTERS = string.ascii_lowercase
UPPERCASE_LETTERS = string.ascii_uppercase
DIGITS = string.digits

skipping_characters = ['(', '=', '<', '>', '?']
PRIVATE_MOTIF_LABEL = '?'
INDEX_TO_CHR = list(LOWERCASE_LETTERS) + list(UPPERCASE_LETTERS) + list(DIGITS)
INDEX_TO_CHR.extend([chr(x) for x in range(33, 127) if chr(x) not in skipping_characters and chr(x) not in INDEX_TO_CHR])

DNA_CHARACTERS = {'A', 'C', 'G', 'T'}


def read_fasta(fasta_file):
    headers = []
    sequences = []
    sequence = ""
    header = ""
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                if len(sequence) > 0:
                    headers.append(header)
                    sequences.append(sequence)
                header = line.strip()[1:]
                sequence = ""
            else:
                sequence += line.strip().upper()
        headers.append(header)
        sequences.append(sequence)

    return headers, sequences


def get_motif_counter(decomposed_vntrs):
    motif_counter = Counter()
    for decomposed_vntr in decomposed_vntrs:
        motif_counter.update(Counter(decomposed_vntr))

    return motif_counter


def is_emitting_state(state_name):
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
    for s in sequence:
        if s not in DNA_CHARACTERS:
            return False
    return True


def sort(aligned_vntrs, sample_ids, method='lexicographically'):
    if method == 'lexicographically':
        return zip(*sorted(list(zip(sample_ids, aligned_vntrs)), key=lambda x: x[0]))
    else:
        raise ValueError("Please check the sorting method. {}".format(method))


def get_levenshtein_distance(s1, s2):
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
                    cost += len(seq2[i])
                else:
                    cost += len(seq1[i])
            # cost += 1
    return cost


def calculate_cost(alinged_vntrs, alphabet_to_motif):
    # same length, start with simple cost
    total_cost = 0
    for i in range(len(alinged_vntrs) - 1):
        total_cost += _calculate_cost(alinged_vntrs[i], alinged_vntrs[i + 1], alphabet_to_motif)

    return total_cost


def sort_with_simulated_annealing(seq_list, alphabet_to_motif, sample_ids):
    initial_cost = calculate_cost(seq_list, alphabet_to_motif)
    T = 10000

    iteration = 0
    while True:
        iteration += 1
        print("T:", T)
        if T <= 1:
            break
        print("iteration", iteration)
        not_changed_count = 0

        for pair in itertools.combinations(range(len(seq_list)), 2):
            first, second = pair

            # only compare the cost before and after changing the order
            current_cost = 0
            after_cost = 0
            # Flanking cost for the first sequence
            current_cost += _calculate_cost(seq_list[first], seq_list[first + 1], alphabet_to_motif)
            after_cost += _calculate_cost(seq_list[second], seq_list[first + 1], alphabet_to_motif)

            current_cost += _calculate_cost(seq_list[first], seq_list[first - 1], alphabet_to_motif)
            after_cost += _calculate_cost(seq_list[second], seq_list[first - 1], alphabet_to_motif)

            # Flanking cost for the second sequence
            current_cost += _calculate_cost(seq_list[second], seq_list[second - 1], alphabet_to_motif)
            after_cost += _calculate_cost(seq_list[first], seq_list[second - 1], alphabet_to_motif)

            if second == len(seq_list) - 1:
                current_cost += _calculate_cost(seq_list[second], seq_list[0], alphabet_to_motif)
                after_cost += _calculate_cost(seq_list[first], seq_list[0], alphabet_to_motif)
            else:
                current_cost += _calculate_cost(seq_list[second], seq_list[second + 1], alphabet_to_motif)
                after_cost += _calculate_cost(seq_list[first], seq_list[second + 1], alphabet_to_motif)

            if after_cost < current_cost:
                seq_list[first], seq_list[second] = seq_list[second], seq_list[first]
                sample_ids[first], sample_ids[second] = sample_ids[second], sample_ids[first]
            else:
                prob = np.exp(-(after_cost - current_cost)*1000/T)
                if prob > np.random.uniform(low=0.0, high=1.0):
                    # swap
                    # print("swap because of prob {}".format(prob), "after cost", after_cost, "cur cost", current_cost)
                    seq_list[first], seq_list[second] = seq_list[second], seq_list[first]
                    sample_ids[first], sample_ids[second] = sample_ids[second], sample_ids[first]
                else:
                    not_changed_count += 1

        T *= 0.90
        # if not_changed_count == len(seq_list):
        #     break

    print("initial cost", initial_cost)
    print("after cost", calculate_cost(seq_list, alphabet_to_motif))
    return seq_list