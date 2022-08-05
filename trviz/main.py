from typing import Dict, List

from trviz.decomposer import TandemRepeatDecomposer
from trviz.motif_aligner import MotifAligner
from trviz.visualization import TandemRepeatVisualizer

import numpy as np

from trviz.utils import sort_lexicographically, read_fasta


class TandemRepeatVizWorker:

    def __init__(self):
        self.decomposer = TandemRepeatDecomposer()
        self.motif_aligner = MotifAligner()
        self.visualizer = TandemRepeatVisualizer()

    def generate_tr_plot(self,
                         tr_sequences,
                         motifs,
                         vid,
                         sample_ids,
                         private_motif_threshold,
                         figure_size=None):

        print("VID: {}".format(vid))
        print("Motifs: {}".format(motifs))

        fasta_out = open("{}_seq.fasta".format(vid), "w")

        decomposed_vntrs = []
        # Test with one motif (for faster decomposition)
        motifs = motifs[0]
        print("Motifs used for decomposition: {}".format(motifs))

        for i, tr_sequence in enumerate(tr_sequences):
            print("Processing: {}, {}".format(i, tr_sequence))
            fasta_out.write(">{}\n{}\n".format(i, tr_sequence))
            decomposed_vntrs.append(self.decomposer.decompose_dp(tr_sequence, motifs, verbose=False))

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
        # TODO lable motifs - should be in the utils?
        labeled_vntrs, motif_to_alphabet, alphabet_to_motif = self.decomposer.label_motifs(decomposed_vntrs,
                                                                                           private_motif_threshold,
                                                                                           auto=True)
        with open("{}_motif_map.txt".format(vid), "w") as f:
            for motif in motif_to_alphabet:
                f.write("{}\t{}\n".format(motif, motif_to_alphabet[motif]))

        print("Motif to alphabet dict", motif_to_alphabet)
        print("Alphabet dict", alphabet_to_motif)
        print("Labeled vntrs", labeled_vntrs)

        # Align VNTRs
        sample_ids, aligned_vntrs = self.motif_aligner.align(sample_ids, labeled_vntrs, vid, tool="mafft")
        print(aligned_vntrs)

        if len(aligned_vntrs) == 0:
            return

        # Sort TR sequences
        sorted_aligned_vntrs, sample_ids = sort_lexicographically(aligned_vntrs, sample_ids)
        print("sorted ids", sample_ids)

        # Visualization
        self.visualizer.plot(sorted_aligned_vntrs,
                        figure_size=figure_size,
                        outfolder="../long_vntr_plots", outname=str(vid) + "_annealing",
                        dpi=300,
                        xticks=sample_ids,
                        xtick_degrees=90,
                        hide_yticks=False)


if __name__ == "__main__":
    # test
    temp_output_name = "../temp/temp_output_ART1.fa"
    aligned_vntrs = []
    sample_ids = []
    tr_seq = None
    with open(temp_output_name, "r") as f:
        for line in f:
            if line.startswith(">"):
                sample_ids.append(line.strip()[1:])
                if tr_seq is not None:
                    aligned_vntrs.append(tr_seq)
                tr_seq = ""
            else:
                tr_seq += line.strip()

    if len(tr_seq) > 0:
        aligned_vntrs.append(tr_seq)

    # Sort TR sequences
    sorted_aligned_vntrs, sample_ids = sort_lexicographically(aligned_vntrs, sample_ids)
    print("sorted ids", sample_ids)

    # Visualization
    visualizer = TandemRepeatVisualizer()
    visualizer.plot(sorted_aligned_vntrs,
                    figure_size=(10, 10),
                    outfolder="../long_vntr_plots", outname=str('ART1') + "_annealing",
                    dpi=50,
                    xticks=sample_ids,
                    xtick_degrees=90,
                    hide_yticks=True,
                    debug=True,
                    )
    exit(1)

    #
    tr_visualizer = TandemRepeatVizWorker()
    from trviz.utils import INDEX_TO_CHR

    # load tr sequences - using fasta (input)
    fasta_file = "../Paul_data/1kg_ART1_alleles.txt"
    hedaers, tr_sequences = read_fasta(fasta_file)
    vid = "ART1"
    motifs = ['CCCATCAGACATCAAGTCTTGTAAATTCCACCTCCTACCTCTGTCTGTCCTCACTGCCATTCT']
    PRIVATE_MOTIF_THRESHOLD = 0

    tr_visualizer.generate_tr_plot(tr_sequences,
                                   motifs,
                                   vid,
                                   hedaers,
                                   PRIVATE_MOTIF_THRESHOLD,
                                   figure_size=(10, 15)
                                   )
    # tr_sequences,
    # motifs,
    # vid,
    # sample_ids,
    # headers,
    # private_motif_threshold,
    # figure_size
