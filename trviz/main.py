from typing import Dict, List
import sys
sys.path.insert(0, './')

from trviz.decomposer import TandemRepeatDecomposer
from trviz.motif_aligner import MotifAligner
from trviz.visualization import TandemRepeatVisualizer
from trviz.utils import write_motif_map
from trviz.utils import sort_lexicographically

import numpy as np


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
        labeled_vntrs, motif_to_alphabet, alphabet_to_motif, motif_counter = self.decomposer.label_motifs(
                                                                                           decomposed_vntrs,
                                                                                           private_motif_threshold,
                                                                                           auto=True)
        motif_map_file = f"{vid}_motif_map.txt"
        write_motif_map(motif_map_file, motif_to_alphabet, motif_counter)

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
                             output_name=f"long_vntr_plots/{str(vid)}",
                             dpi=300,
                             sample_ids=sample_ids,
                             xtick_degrees=90,
                             hide_yticks=False)
