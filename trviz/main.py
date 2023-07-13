import sys
import os
from typing import List, Tuple, Dict

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
                        tr_id: str,
                        sample_ids: List[str],
                        tr_sequences: List[str],
                        motifs: List[str],

                        skip_alignment: bool = False,
                        rearrangement_method: str = 'clustering',
                        sample_order_file: str = None,
                        output_dir: str = "./",

                        # Figure parameters
                        figure_size: Tuple[int, int] = None,
                        output_name: str = None,
                        dpi: int = 300,
                        hide_xticks: bool = False,
                        hide_yticks: bool = False,
                        hide_dendrogram: bool = True,
                        population_data: str = None,
                        allele_as_row: bool = True,
                        xlabel_size: int = 8,
                        ylabel_size: int = 8,
                        xlabel_rotation: int = 0,
                        ylabel_rotation: int = 0,
                        private_motif_color: str = 'black',
                        frame_on: Dict[str, bool] = None,
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

        :param tr_id: a ID for the tandem repeat
        :param sample_ids: a list of sample IDs corresponding to the tandem repeat sequences
        :param tr_sequences: a list of tandem repeat sequences
        :param motifs: a list of motifs to be used for decomposition

        :param skip_alignment: if true, skip the multiple sequence alignment
        :param rearrangement_method: options: {'clustering' (default), 'name', 'motif_count',
                                               'simulated_annealing', 'manually'}
        :param sample_order_file: a file containing sample order
        :param output_dir: base directory for output files

        :param figure_size: figure size
        :param output_name: output file name
        :param dpi: DPI for the plot
        :param hide_xticks: if true, hide xticks
        :param hide_yticks: if true, hide yticks
        :param hide_dendrogram: if true, hide the dendrogram
        :param population_data: population data file name
        :param allele_as_row: if true, plot allele as row (default is true)
        :param xlabel_size: x label size (default is 8)
        :param ylabel_size: y label size (default is 8)
        :param xlabel_rotation: x label rotation (default is 0)
        :param ylabel_rotation: y label rotation (default is 0)
        :param private_motif_color: the color for private motifs. Default is black
        :param frame_on: a dictionary mapping sample to frame on.
                         Default is {'top': False, 'bottom': True, 'right': False, 'left': True}
        :param verbose: if true, output detailed information
        """

        if len(sample_ids) != len(tr_sequences):
            raise ValueError("The number of sample IDs and the number of sequences are different.")

        if rearrangement_method == 'manually':
            if sample_order_file is None:
                raise ValueError("Please specify the sample order file for manual rearrangement.")
            if not os.path.exists(sample_order_file):
                raise ValueError("The specified sample order file does not exist.")

        if verbose:
            print("VID: {}".format(tr_id))
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
                                                motif_map_file=f"{output_dir}/{tr_id}_motif_map.txt",
                                                auto=True)

        # 3. Align motifs
        if skip_alignment:
            if verbose:
                print("Skip alignment step")
            aligned_trs = add_padding(encoded_trs)
            sorted_sample_ids = sample_ids
        else:
            if verbose:
                print("Alignment")
            score_matrix = get_score_matrix(self.motif_encoder.symbol_to_motif)
            sorted_sample_ids, aligned_trs = self.motif_aligner.align(sample_ids,
                                                                      encoded_trs,
                                                                      tr_id,
                                                                      score_matrix=score_matrix,
                                                                      output_dir=output_dir)

        # 4. Re-arrangement
        if rearrangement_method is not None and rearrangement_method != 'clustering':
            sorted_sample_ids, aligned_trs = sort(aligned_trs,
                                                  sorted_sample_ids,
                                                  self.motif_encoder.symbol_to_motif,
                                                  sample_order_file,
                                                  rearrangement_method)

        # 5. Visualization
        if verbose:
            print("Visualization")
        self.visualizer.trplot(aligned_labeled_repeats=aligned_trs,
                               sample_ids=sorted_sample_ids,
                               figure_size=figure_size,
                               output_name=f"{output_dir}/{str(tr_id)}.pdf" if output_name is None else output_name,
                               dpi=dpi,
                               sort_by_clustering=True if rearrangement_method == 'clustering' else False,
                               hide_xticks=hide_xticks,
                               hide_yticks=hide_yticks,
                               allele_as_row=allele_as_row,
                               xlabel_size=xlabel_size,
                               ylabel_size=ylabel_size,
                               hide_dendrogram=hide_dendrogram,
                               symbol_to_motif=self.motif_encoder.symbol_to_motif,
                               xlabel_rotation=xlabel_rotation,
                               ylabel_rotation=ylabel_rotation,
                               population_data=population_data,
                               private_motif_color=private_motif_color,
                               frame_on=frame_on
                               )

        # 6. Motif map
        self.visualizer.plot_motif_color_map(self.motif_encoder.symbol_to_motif,
                                             self.motif_encoder.motif_counter,
                                             self.visualizer.symbol_to_color,
                                             file_name=f"{output_dir}/{str(tr_id)}_motif_map.png",
                                             )
