__author__ = "Jonghun Park"
__date__ = 1 / 12 / 22
__email__ = "jonghunpark@ucsd.edu"

import glob

from decomposer import decompose_dp, label_motifs
from motif_aligner import align_motifs
from visualization import trplot

import itertools

from scipy.cluster import hierarchy
# from scipy.spatial.distance import pdist
# from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import numpy as np


def levenshteinDistance(s1, s2):
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
                cost += levenshteinDistance(s1, s2)
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
        total_cost += _calculate_cost(alinged_vntrs[i], alinged_vntrs[i+1], alphabet_to_motif)

    return total_cost


def sort_with_simulated_annealing(seq_list, alphabet_to_motif):
    initial_cost = calculate_cost(seq_list, alphabet_to_motif)
    T = 10000

    iteration = 0
    while True:
        iteration += 1
        T *= 0.90
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

            # if first == 0:
            #     current_cost += _calculate_cost(seq_list[first], seq_list[first - 1])
            #     after_cost += _calculate_cost(seq_list[second], seq_list[first - 1])
            # else:
            #     current_cost += _calculate_cost(seq_list[first], seq_list[first - 1])
            #     after_cost += _calculate_cost(seq_list[second], seq_list[first - 1])

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
                # swap
                print("SWAPED", pair, after_cost, current_cost)
                if first != 0:
                    print(seq_list[first-1])
                print(seq_list[first])
                print(seq_list[first+1])

                print(seq_list[second-1])
                print(seq_list[second])
                if second != len(seq_list) - 1:
                    print(seq_list[second+1])
                seq_list[first], seq_list[second] = seq_list[second], seq_list[first]
            else:
                prob = np.exp(-(after_cost - current_cost)*1000/T)
                if prob > np.random.uniform(low=0.0, high=1.0):
                    # swap
                    # print("swap because of prob {}".format(prob), "after cost", after_cost, "cur cost", current_cost)
                    seq_list[first], seq_list[second] = seq_list[second], seq_list[first]
                else:
                    not_changed_count += 1

        # if not_changed_count == len(seq_list):
        #     break

    print("initial cost", initial_cost)
    print("after cost", calculate_cost(seq_list, alphabet_to_motif))
    return seq_list


def sort_lexicographically(aligned_vntrs, sample_ids):
    return zip(*sorted(list(zip(aligned_vntrs, sample_ids)), key=lambda x: x[0]))


def is_polymorphic(tr_seqs, threshold=0.1):
    return True if len(set(tr_seqs)) / float(len(tr_seqs)) > threshold else False


def process_vntr(tr_sequences, motifs, vid, sample_ids, private_motif_threshold):
    print("VID: {}".format(vid))
    print("Motifs: {}".format(motifs))

    ### If already have the alignment output
    # temp_output_name = "temp/temp_output_{}.fa".format(vid)
    # # TEST
    # aligned_vntrs = []
    # tr_seq = None
    # with open(temp_output_name, "r") as f:
    #     for line in f:
    #         if line.startswith(">"):
    #             if tr_seq is not None:
    #                 aligned_vntrs.append(tr_seq)
    #             tr_seq = ""
    #         else:
    #             tr_seq += line.strip()
    #
    # if len(tr_seq) > 0:
    #     aligned_vntrs.append(tr_seq)
    #
    # # Sort TR sequences
    # sorted_aligned_vntrs, sample_ids = sort_lexicographically(aligned_vntrs, sample_ids)
    #
    # # Visualization
    # trplot(sorted_aligned_vntrs,
    #        outfolder="long_vntr_plots", outname=str(vid) + "_annealing",
    #        dpi=100,
    #        xticks=sample_ids,
    #        xtick_degrees=90,
    #        hide_yticks=False)
    # return

    # Test for clustering

    # 1. Test edit_distance
    # Calculate pairwise edit-distance
    # pairwise_distances = []
    # for i in range(len(tr_sequences) - 1):
    #     for j in range(i + 1, len(tr_sequences)):
    #         dist = levenshteinDistance(tr_sequences[i], tr_sequences[j])
    #         pairwise_distances.append(dist)
    #
    # print(pairwise_distances)
    # ytdist = np.array(pairwise_distances)
    # Z = hierarchy.linkage(ytdist, 'single')
    # plt.figure()
    # dn = hierarchy.dendrogram(Z)
    # plt.show()
    # exit(1)

    # TEST FASTA
    fasta_out = open("{}_seq.fasta".format(vid), "w")

    decomposed_vntrs = []
    # Test with one motif (for faster decomposition)
    motifs = motifs[0]
    print("Motifs used for decomposition: {}".format(motifs))

    for i, tr_sequence in enumerate(tr_sequences):
        print("Processing: {}, {}".format(i, tr_sequence))
        fasta_out.write(">{}\n{}\n".format(i, tr_sequence))
        decomposed_vntrs.append(decompose_dp(tr_sequence, motifs, verbose=False))
    print("decomposed VNTRs:", decomposed_vntrs)

    # Sort by sequence length
    seq_array = np.array(tr_sequences)
    length_array = np.array([len(s) for s in tr_sequences])
    length_indices = np.argsort(length_array)

    # Sort by motif length
    motif_length_array = np.array([len(s) for s in decomposed_vntrs])
    motif_length_indices = np.argsort(motif_length_array)

    fasta_out.close()

    # Label motifs
    labeled_vntrs, motif_to_alphabet, alphabet_to_motif = label_motifs(decomposed_vntrs, private_motif_threshold, auto=True)
    with open("{}_motif_map.tsv", "w") as f:
        for motif in motif_to_alphabet:
            f.write("{}\t{}\n".format(motif, motif_to_alphabet[motif]))

    print("Motif to alphabet dict", motif_to_alphabet)
    print("Alphabet dict", alphabet_to_motif)
    print("Labeled vntrs", labeled_vntrs)

    # Align VNTRs
    aligned_vntrs = align_motifs(labeled_vntrs, vid, method="mafft")
    print(aligned_vntrs)

    if len(aligned_vntrs) == 0:
        return

    # Sort TR sequences
    sorted_aligned_vntrs, sample_ids = sort_lexicographically(aligned_vntrs, sample_ids)

    # Visualization
    trplot(sorted_aligned_vntrs,
           outfolder="long_vntr_plots", outname=str(vid) + "_annealing",
           dpi=100,
           xticks=sample_ids,
           xtick_degrees=90,
           hide_yticks=False)
    return

    # IMPORTANT!!!!!!!!!!!!!!
    # THIS ALL DOESN'T WORK. MUSCLE CHANGES THE ORDER OR SEQ.
    # IMPORTANT!!!!!!!!!!!!!!

    # TEST parsing guide tree
    # index_order = []
    # with open("690585_seq.guidetree", "r") as f:
    #     for line in f:
    #         if ":" in line and not line.startswith(")"):
    #             index = int(line.split(":")[0])
    #             index_order.append(index)
    #
    # assert len(index_order) == len(aligned_vntrs)

    # Sort by guide tree
    # sorted_aligned_vntrs = []
    # for i in index_order:
    #     sorted_aligned_vntrs.append(aligned_vntrs[i])
    # aligned_vntrs = sorted_aligned_vntrs

    # Sort by tr seuqnece length
    # sorted_aligned_vntrs = []
    # for i in length_indices:
    #     sorted_aligned_vntrs.append(aligned_vntrs[i])
    # aligned_vntrs = sorted_aligned_vntrs

    # Sort by motif length
    motif_len = []
    for vntr in aligned_vntrs:
        length = 0
        for v in vntr:
            if v != "-":
                length += 1
        motif_len.append(length)

    motif_length_array = np.array(motif_len)
    motif_length_indices = np.argsort(motif_length_array)

    # Sort by label length
    sorted_aligned_vntrs = []
    for i in motif_length_indices:
        sorted_aligned_vntrs.append(aligned_vntrs[i])
    aligned_vntrs = sorted_aligned_vntrs

    # Get length boundaries
    motif_len = []
    for vntr in aligned_vntrs:
        length = 0
        for v in vntr:
            if v != "-":
                length += 1
        motif_len.append(length)

    prev_length = 0
    length_boundaries = []
    for i, vntr in enumerate(aligned_vntrs):
        length = 0
        for v in vntr:
            if v != "-":
                length += 1

        if length != prev_length:
            length_boundaries.append(i)
            prev_length = length

    length_boundaries.append(len(aligned_vntrs))  # last one

    # for each length
    sorted_aligned_vntrs = []
    for i in range(len(length_boundaries)-1):
        same_length_motifs = aligned_vntrs[length_boundaries[i]: length_boundaries[i+1]]
        # print("before")
        # print(same_length_motifs)
        # sorted_vntrs = sort_with_simulated_annealing(same_length_motifs, alphabet_to_motif)
        # print("after")
        # print(sorted_vntrs)
        sorted_aligned_vntrs.extend(same_length_motifs)
        # sorted_aligned_vntrs.extend(sorted_vntrs)

    # Simulated annealing
    # print("Simulated annealing")
    # alinged_vntrs = sort_with_simulated_annealing(aligned_vntrs)

    # Visualization
    trplot(sorted_aligned_vntrs,
           outfolder="long_vntr_plots", outname=str(vid) + "_annealing",
           dpi=100,
           alpha=0.5,
           xticks=sample_ids,
           xtick_degrees=90,
           hide_yticks=False)


if __name__ == "__main__":
    # load vntr data


    # target_vntr_files = glob.glob("./tr_sequences/*sequences.tsv")
    # Error VNTR
    target_vntr_files = glob.glob("tr_sequences/385941_tr_sequences.tsv")

    # target_vntr_files = glob.glob("./tr_sequences/830537_tr_sequences.tsv")
    # target_vntr_files = glob.glob("./tr_sequences/329655_tr_sequences_edited.tsv")
    # target_vntr_files = glob.glob("./tr_sequences/329655_tr_sequences_selected.tsv")

    target_vntr_files = glob.glob("course_tr_sequences/*.tsv")
    # target_vntr_files = glob.glob("tr_sequences/329655_tr_sequences.tsv")  # VPS53 41 bp motif
    # target_vntr_files = glob.glob("tr_sequences/830537_tr_sequences.tsv")  # WDR60 38 bp motif

    PRIVATE_MOTIF_THRESHOLD = 0
    print(len(target_vntr_files))

    for target_vntr_file in target_vntr_files:
        print(target_vntr_file)
        vntr_id = int(target_vntr_file.split("/")[-1].split("_")[0])
        motifs = []
        tr_seqs = []
        haplotype_ids = []
        with open(target_vntr_file , "r") as f:
            for line in f:
                if line.startswith("#"):
                    if "Motifs" in line:
                        motifs = line.strip().split(" ")[-1].split("\t")
                else:
                    split = line.strip().split("\t")
                    if len(split) == 1:  # no tr seq
                        continue
                    filename, tr_seq = split
                    haplotype_ids.append(filename.split("/")[-1].split("_")[1])
                    tr_seqs.append(tr_seq)

        print("number of haplotypes", len(haplotype_ids))
        print("number of trseqs", len(tr_seqs))

        if not is_polymorphic(tr_seqs):
            print("Not polymorphic")
            continue

        if len(tr_seqs) < 20:
            print("Too small number of TR sequences")
            continue

        process_vntr(tr_seqs, motifs, vntr_id, haplotype_ids, private_motif_threshold=PRIVATE_MOTIF_THRESHOLD)

    exit(1)


    # For each VNTR
    trseq_input_file = "tr_sequences.txt"
    MAX_SAMPLE_COUNT = 1000
    with open(trseq_input_file, "r") as f:
        vid = None
        for line in f:
            if "VID: " in line:
                # Process and visualize the VNTR
                if vid is not None:
                    if vid == 353349:  # This has more than 21 motifs.
                    # if vid == 690585:  # test the weird case first
                        process_vntr(tr_seqs, motifs, vid, sample_ids)
                        exit(1)
                    # process_vntr(tr_seqs, motifs, vid, sample_ids)

                # Initialize variables
                vid = int(line.strip().split(" ")[-1])
                tr_seqs = []
                sample_ids = []
                sample_count = 0
            elif "Genotype: " in line:
                genotype = line.strip().split(" ")[-1]
            elif "Motifs: " in line:
                motifs = line.strip().split(" ")[-1].split("\t")
                print("Motifs:", motifs)
            elif "Sample: " in line:
                sample = line.strip().split(" ")[-1]
                sample_count += 1
            else:  # Format: Copy\tTandemRepeatSequence
                try:
                    copy_number, tr_seq = line.strip().split("\t")
                except:
                    print("ERROR found in line: ", line)
                    exit(1)
                if sample_count < MAX_SAMPLE_COUNT:
                    sample_ids.append(sample)  # One sample can have up to 2 alleles
                    tr_seqs.append(tr_seq)

    # Process the last VNTR
    if vid is not None:
        process_vntr(tr_seqs, motifs, vid, sample_ids)
