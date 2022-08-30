import sys
sys.path.insert(0, './')

from trviz.decomposer import TandemRepeatDecomposer
from trviz.motif_encoder import MotifEncoder
from trviz.motif_aligner import MotifAligner
from trviz.visualizer import TandemRepeatVisualizer
from trviz.utils import sort


class TandemRepeatVizWorker:

    def __init__(self):
        self.decomposer = TandemRepeatDecomposer()
        self.motif_encoder = MotifEncoder()
        self.motif_aligner = MotifAligner()
        self.visualizer = TandemRepeatVisualizer()

    def generate_tr_plot(self,
                         vntr_id,
                         sample_ids,
                         tr_sequences,
                         motifs,
                         figure_size=None,
                         rearragement=None,
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
                print(f"Decomposing TR sequence {i}")
            decomposed_vntrs.append(self.decomposer.decompose(tr_sequence, motifs))

        # 2. Encoding
        encoded_vntrs = self.motif_encoder.encode(decomposed_vntrs,
                                                  motif_map_file=f"{vntr_id}_motif_map.txt",
                                                  auto=True)

        # 3. Align motifs
        sample_ids, aligned_vntrs = self.motif_aligner.align(sample_ids, encoded_vntrs, vntr_id)

        # 4. Sorting
        if rearragement is not None:
            sample_ids, aligned_vntrs = sort(aligned_vntrs, sample_ids)

        # 5. Visualization
        self.visualizer.plot(aligned_vntrs,
                             figure_size=figure_size,
                             output_name=f"long_vntr_plots/{str(vntr_id)}",
                             dpi=300,
                             sample_ids=sample_ids,
                             xtick_degrees=90,
                             hide_yticks=False)
