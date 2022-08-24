import sys
sys.path.insert(0, './')

from trviz.decomposer import TandemRepeatDecomposer
from trviz.motif_aligner import MotifAligner
from trviz.visualization import TandemRepeatVisualizer
from trviz.utils import write_motif_map
from trviz.utils import sort_lexicographically


class TandemRepeatVizWorker:

    def __init__(self):
        self.decomposer = TandemRepeatDecomposer()
        self.motif_aligner = MotifAligner()
        self.visualizer = TandemRepeatVisualizer()

    def generate_tr_plot(self,
                         tr_sequences,
                         motifs,
                         vntr_id,
                         sample_ids,
                         private_motif_threshold=0,
                         figure_size=None,
                         verbose=True,
                         ):

        if verbose:
            print("VID: {}".format(vntr_id))
            print("Motifs: {}".format(motifs))
            print(f"{len(tr_sequences)} sequences")

        # 1. Decomposition
        decomposed_vntrs = []
        for i, tr_sequence in enumerate(tr_sequences):
            if verbose:
                print(f"Decomposing TR {i}")
            decomposed_vntrs.append(self.decomposer.decompose_dp(tr_sequence, motifs, verbose=False))
        if verbose:
            print(f"Decomposed VNTRs: {decomposed_vntrs}")

        # 2. Labeling
        labeled_vntrs, motif_to_alphabet, alphabet_to_motif, motif_counter = self.decomposer.label_motifs(
                                                                                           decomposed_vntrs,
                                                                                           private_motif_threshold,
                                                                                           auto=True)
        # Write motif map
        motif_map_file = f"{vntr_id}_motif_map.txt"
        write_motif_map(motif_map_file, motif_to_alphabet, motif_counter)

        if verbose:
            print("Labeled vntrs", labeled_vntrs)

        # 3. Align motifs
        sample_ids, aligned_vntrs = self.motif_aligner.align(sample_ids, labeled_vntrs, vntr_id, tool="mafft")

        # Sort TR sequences
        sorted_aligned_vntrs, sample_ids = sort_lexicographically(aligned_vntrs, sample_ids)

        # Visualization
        self.visualizer.plot(sorted_aligned_vntrs,
                             figure_size=figure_size,
                             output_name=f"long_vntr_plots/{str(vntr_id)}",
                             dpi=300,
                             sample_ids=sample_ids,
                             xtick_degrees=90,
                             hide_yticks=False)
