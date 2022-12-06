import sys
from typing import List, Tuple

sys.path.insert(0, './')

from trviz.decomposer import Decomposer
from trviz.motif_encoder import MotifEncoder
from trviz.motif_aligner import MotifAligner
from trviz.visualizer import TandemRepeatVisualizer
from trviz.utils import sort
from trviz.utils import add_padding
from trviz.utils import get_score_matrix


class TandemRepeatVizWorker:

    def __init__(self):
        self.decomposer = Decomposer()
        self.motif_encoder = MotifEncoder()
        self.motif_aligner = MotifAligner()
        self.visualizer = TandemRepeatVisualizer()

    def generate_trplot(self,
                        vntr_id: str,
                        sample_ids: List[str],
                        tr_sequences: List[str],
                        motifs: List[str],
                        figure_size: Tuple[int, int] = None,
                        rearrangement_method: str = 'clustering',
                        hide_dendrogram: bool = False,
                        skip_alignment: bool = False,
                        output_dir: str = "./",
                        verbose: bool = True,
                        ):
        """
        A method to generate a plot of tandem repeat motif composition.
        This executes the following modules sequentially and finally output the plot.
        1. Decomposition
        2. Encoding
        3. Alignment
        4. Rearrangement (if specified)
        5. Visualization
        For detail, please check out each module

        :param vntr_id: a ID for the tandem repeat
        :param sample_ids: a list of sample IDs corresponding to the tandem repeat sequences
        :param tr_sequences: a list of tandem repeat sequences
        :param motifs: a list of motifs to be used for decomposition
        :param figure_size: figure size
        :param rearrangement_method: options: {'clustering' (default), 'lexicographically', 'simulated_annealing'}
        :param hide_dendrogram: if True, hide dendrogram
        :param skip_alignment: if true, skip the multiple sequence alignment
        :param output_dir: base directory for output files
        :param verbose: if true, output detailed information
        """

        if verbose:
            print("VID: {}".format(vntr_id))
            print("Motifs: {}".format(motifs))
            print(f"Loaded {len(tr_sequences)} tandem repeat sequences")
            print("Decomposing TR sequences")

        # 1. Decomposition
        decomposed_trs = []
        for i, tr_sequence in enumerate(tr_sequences):
            if verbose:
                from trviz.utils import print_progress_bar
                print_progress_bar(i + 1, len(tr_sequences))
            decomposed_trs.append(self.decomposer.decompose(tr_sequence, motifs))

        # 2. Encoding
        if verbose:
            print("Encoding")
        encoded_trs = self.motif_encoder.encode(decomposed_trs,
                                                motif_map_file=f"{output_dir}/{vntr_id}_motif_map.txt",
                                                auto=True)

        # 3. Align motifs
        if skip_alignment:
            if verbose:
                print("Skip alignment step")
            aligned_trs = add_padding(encoded_trs)
        else:
            if verbose:
                print("Alignment")
            score_matrix = get_score_matrix(self.motif_encoder.symbol_to_motif)
            sorted_sample_ids, aligned_trs = self.motif_aligner.align(sample_ids,
                                                                      encoded_trs,
                                                                      vntr_id,
                                                                      score_matrix=score_matrix,
                                                                      output_dir=output_dir)

        # 4. Re-arrangement
        if rearrangement_method is not None and rearrangement_method != 'clustering':
            sorted_sample_ids, aligned_trs = sort(aligned_trs,
                                                  sorted_sample_ids,
                                                  self.motif_encoder.symbol_to_motif,
                                                  rearrangement_method)

        # 5. Visualization
        if verbose:
            print("Visualization")
        self.visualizer.trplot(aligned_labeled_repeats=aligned_trs,
                               sample_ids=sorted_sample_ids,
                               figure_size=figure_size,
                               output_name=f"{output_dir}/{str(vntr_id)}",
                               dpi=500,
                               xtick_degrees=90,
                               sort_by_clustering=True if rearrangement_method == 'clustering' else False,
                               hide_dendrogram=hide_dendrogram,
                               symbol_to_motif=self.motif_encoder.symbol_to_motif,
                               )

        # 6. Motif map
        self.visualizer.plot_motif_color_map(self.motif_encoder.symbol_to_motif,
                                             self.motif_encoder.motif_counter,
                                             self.visualizer.symbol_to_color,
                                             f"{output_dir}/{str(vntr_id)}_motif_map.png",
                                             )
