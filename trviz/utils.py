import string

LOWERCASE_LETTERS = string.ascii_lowercase
UPPERCASE_LETTERS = string.ascii_uppercase
DIGITS = string.digits

skipping_characters = ['(', '=', '<', '>']
INDEX_TO_CHR = list(LOWERCASE_LETTERS) + list(UPPERCASE_LETTERS) + list(DIGITS)
INDEX_TO_CHR.extend([chr(x) for x in range(33, 127) if chr(x) not in skipping_characters and chr(x) not in INDEX_TO_CHR])

DNA_CHARACTERS = {'A', 'C', 'G', 'T'}


def read_fasta(fasta_file):
    headers = []
    sequences = []
    sequence = ""
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


def sort_lexicographically(aligned_vntrs, sample_ids):
    return zip(*sorted(list(zip(aligned_vntrs, sample_ids)), key=lambda x: x[0]))
