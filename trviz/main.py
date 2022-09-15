import sys
from typing import List, Tuple

sys.path.insert(0, './')

from trviz.decomposer import Decomposer
from trviz.motif_encoder import MotifEncoder
from trviz.motif_aligner import MotifAligner
from trviz.visualizer import TandemRepeatVisualizer
from trviz.utils import sort


class TandemRepeatVizWorker:

    def __init__(self):
        self.decomposer = Decomposer()
        self.motif_encoder = MotifEncoder()
        self.motif_aligner = MotifAligner()
        self.visualizer = TandemRepeatVisualizer()

    def generate_tr_plot(self,
                         vntr_id: str,
                         sample_ids: List[str],
                         tr_sequences: List[str],
                         motifs: List[str],
                         figure_size: Tuple[int, int] = None,
                         rearragement: str = None,
                         verbose: bool = True,
                         ):
        """
        A method to generate TandemRepeat plot.
        This executes the following modules sequentially and finally output the plot.
        1. Decomposition
        2. Encoding
        3. Alignment
        4. Sorting (if specified)
        5. Visualization
        For detail, please check out each module

        :param vntr_id: a ID for the tandem repeat
        :param sample_ids: a list of sample IDs corresponding to the tandem repeat sequences
        :param tr_sequences: a list of tandem repeat sequences
        :param motifs: a list of motifs to be used for decomposition
        :param figure_size: figure size
        :param rearragement: re-arrangement method (default is sorting by lexicographically)
        :param verbose: if true, output detailed information
        """

        if verbose:
            print("VID: {}".format(vntr_id))
            print("Motifs: {}".format(motifs))
            print(f"Loaded {len(tr_sequences)} tandem repeat sequences")

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
                             output_name=f"{str(vntr_id)}",
                             dpi=300,
                             sample_ids=sample_ids,
                             xtick_degrees=90,
                             hide_yticks=False)
