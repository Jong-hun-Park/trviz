from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from trviz.utils import PRIVATE_MOTIF_LABEL


class TandemRepeatVisualizer:

    def __init__(self):
        self.symbol_to_color = None

    @staticmethod
    def encode_tr_sequence(labeled_motifs):
        decomposed_motifs = []
        for motif in labeled_motifs:
            encoded = []
            for m in motif:
                if m == '-':
                    encoded.append(0)
                else:
                    encoded.append(ord(m) - 64)  # 1-base start
            decomposed_motifs.append(encoded)

        return decomposed_motifs

    @staticmethod
    def _get_unique_labels(aligned_repeats):
        unique_repeats = set()
        for rs in aligned_repeats:
            for r in rs:
                unique_repeats.add(r)

        return sorted(unique_repeats)

    @staticmethod
    def plot_motif_color_map(symbol_to_motif, motif_counter, symbol_to_color, file_name,
                             figure_size=None,
                             box_size=1,
                             box_margin=0.1,
                             label_size=None,
                             dpi=300):
        """ Generate a plot for motif color map """

        if box_margin < 0 or box_margin > box_size:
            raise ValueError(f"Box margin should be between 0 and {box_size}.")

        if symbol_to_color is None:
            raise ValueError("Symbol to color is not set.")

        if figure_size is not None:
            fig, ax = plt.subplots(figsize=figure_size)
        else:
            max_motif_length = len(max(symbol_to_motif.values(), key=len))
            w = len(symbol_to_motif) // 10 + 5 if len(symbol_to_motif) > 50 else len(symbol_to_motif) // 5 + 1
            h = max_motif_length // 4 + 1 if max_motif_length > 50 else max_motif_length // 3 + 1
            if h * dpi > 2 ** 16:
                h = 2 ** 15 // dpi * 0.75  # "Weight and Height must be less than 2^16"
            fig, ax = plt.subplots(figsize=(w, h))
        ax.set_aspect('equal')

        max_motif_length = len(max(symbol_to_motif.values(), key=len))
        # TODO: sort by the frequency
        # TODO: figure and font size
        yticklabels = []
        y_base_position = 0
        has_gap_in_symbol_to_color = False
        for (symbol, color) in symbol_to_color.items():
            if symbol == '-':
                has_gap_in_symbol_to_color = True
                continue
            box_position = [box_margin, y_base_position + box_margin]
            box_width = box_size - 2 * box_margin
            box_height = box_size - 2 * box_margin
            ax.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                       linewidth=0,
                                       facecolor=color,
                                       edgecolor="white"))
            y_base_position += 1

            # text_position = [box_position[0] + 2.5 * box_width / 2, box_position[1] + box_height / 2 - 0.1]
            motif = symbol_to_motif[symbol]
            text = f"{motif:<{max_motif_length + 2}}" + f"{motif_counter[motif]:>5}"
            # ax.text(x=text_position[0], y=text_position[1], s=text, fontname="monospace")
            yticklabels.append(text)

        ax.yaxis.tick_right()
        ax.set_ylim(top=len(symbol_to_color) - 1)
        yticks_count = len(symbol_to_color)
        if has_gap_in_symbol_to_color:
            yticks_count = len(symbol_to_color) - 1
        ax.set_yticks([y * box_size + box_size*0.5 for y in range(yticks_count)])  # minus 1 for gap, '-'
        ax.set_yticklabels(yticklabels)

        # font size
        if label_size is not None:
            ax.tick_params(axis='both', which='major', labelsize=label_size)
            ax.tick_params(axis='both', which='minor', labelsize=label_size)
        # font
        for tick in ax.get_yticklabels():
            tick.set_fontproperties("monospace")

        # Remove yticks
        ax.tick_params(
            axis='y',  # changes apply to the y-axis
            which='both',  # both major and minor ticks are affected
            length=0)  # labels along the bottom edge are on

        # Hide xticks
        ax.set_xticks([])

        # No frame
        ax.set_frame_on(False)

        fig.tight_layout()
        fig.savefig(file_name, dpi=dpi, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

    def trplot(self,
               aligned_labeled_repeats: List[str],
               sample_ids: List[str],
               figure_size: Tuple[int, int] = None,
               output_name: str = None,
               dpi: int = 300,
               alpha: float = 0.6,
               box_line_width: float = 0,
               hide_xticks: bool = False,
               hide_yticks: bool = False,
               symbol_to_motif: Dict[str, str] = None,
               sort_by_clustering: bool = True,
               hide_dendrogram: bool = True,
               population_data: str = None,
               motif_marks: Dict[str, str] = None,
               allele_as_row: bool = True,
               xlabel_size: int = 8,
               ylabel_size: int = 8,
               xlabel_rotation: int = 0,
               ylabel_rotation: int = 0,
               private_motif_color: str = 'black',
               frame_on: Dict[str, bool] = None,
               debug: bool = False
               ):
        """
        Generate a plot showing the variations in tandem repeat sequences.
        A distinct color is assigned to each motif.
        For private motifs (with low frequency), the same color (in black by default) may be assigned.

        :param aligned_labeled_repeats: aligned and encoded tandem repeat sequences.
        :param sample_ids: sample IDs
        :param figure_size: figure size
        :param output_name: output file name
        :param dpi: DPI for the plot
        :param alpha: alpha value for the plot
        :param box_line_width: line width for box edges
        :param hide_xticks: if true, hide xticks
        :param hide_yticks: if true, hide yticks
        :param symbol_to_motif: a dictionary mapping symbols to motifs
        :param sort_by_clustering: if true, sort the samples by clustering
        :param hide_dendrogram: if true, hide the dendrogram
        :param population_data: population data file name
        :param motif_marks: a dictionary mapping sample to motif marks
        :param allele_as_row: if true, plot allele as row (default is true)
        :param xlabel_size: x label size (default is 8)
        :param ylabel_size: y label size (default is 8)
        :param xlabel_rotation: x label rotation (default is 0)
        :param ylabel_rotation: y label rotation (default is 0)
        :param private_motif_color: the color for private motifs. Default is black
        :param frame_on: a dictionary mapping sample to frame on.
                         Default is {'top': False, 'bottom': True, 'right': False, 'left': True}
        :param debug: if true, print verbose information.
        """

        if motif_marks is not None:
            if len(motif_marks) != len(sample_ids):
                raise ValueError("The length of motif_marks must be the same as the number of samples")

        fig, ax_main = plt.subplots(figsize=figure_size)  # width and height in inch
        max_repeat_count = len(aligned_labeled_repeats[0])
        if figure_size is None:
            if allele_as_row:
                h = len(sample_ids) // 5 + 2 if len(sample_ids) > 50 else 5
                w = max_repeat_count // 5 + 2 if max_repeat_count > 50 else max_repeat_count // 5 + 2
            else:
                w = len(sample_ids) // 5 + 2 if len(sample_ids) > 50 else 5
                h = max_repeat_count // 5 + 2 if max_repeat_count > 50 else max_repeat_count // 5 + 2
            if h * dpi > 2**16:
                h = 2**15 // dpi * 0.75  # "Weight and Height must be less than 2^16"
            fig, ax_main = plt.subplots(figsize=(w, h))

        # Add clustering
        if sort_by_clustering:
            if symbol_to_motif is None:
                raise ValueError("symbol_to_motif must be provided when sort_by_clustering is True")
            sorted_sample_ids, sorted_aligned_labeled_repeats = self.add_dendrogram(fig,
                                                                                    aligned_labeled_repeats,
                                                                                    sample_ids,
                                                                                    symbol_to_motif,
                                                                                    hide_dendrogram)
            if debug:
                print("Sort by clustering")
                print('\n'.join(sorted_sample_ids))
        else:
            # Already sorted in some way
            sorted_sample_ids, sorted_aligned_labeled_repeats = sample_ids, aligned_labeled_repeats
            if debug:
                print("No clustering")

        # Set symbol to color map
        unique_labels = self._get_unique_labels(aligned_labeled_repeats)
        unique_label_count = len(unique_labels)
        symbol_to_color = self.get_symbol_to_color_map(alpha, unique_label_count, unique_labels)
        self.set_symbol_to_color_map(symbol_to_color)

        box_height = 1.0
        box_width = 1.0
        for allele_index, allele in enumerate(sorted_aligned_labeled_repeats):
            motif_index = 0
            if allele_as_row:
                box_position = [allele_index, len(sorted_aligned_labeled_repeats) - allele_index - 1]
            else:  # each column is an allele
                box_position = [allele_index, 0]

            for box_index, symbol in enumerate(allele):
                if allele_as_row:
                    box_position[0] = box_width * box_index  # move x position
                else:
                    box_position[1] = box_height * box_index
                fcolor = symbol_to_color[symbol]
                hatch_pattern = None

                if symbol == '-':  # For gaps, color them as white blocks
                    fcolor = (1, 1, 1, 1)
                else:  # Not a gap
                    if motif_marks is not None:
                        motif_mark = motif_marks[sorted_sample_ids[allele_index]][motif_index]
                        if motif_mark == 'I':  # introns
                            hatch_pattern = 'xxx'
                    motif_index += 1

                if symbol == PRIVATE_MOTIF_LABEL:
                    fcolor = private_motif_color
                ax_main.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                                linewidth=box_line_width + 0.1,
                                                facecolor=fcolor,
                                                edgecolor="white",
                                                hatch=hatch_pattern, ))

        # population_data = "../../HPRC_metadata/sample_metadata.tsv"
        if population_data is not None:
            # TODO: need to be fixed when allele_as_row is True (the below works only for allele as column)
            ax_bottom = fig.add_axes([0.125, 0.07, 0.775, 0.03])  # xmin, ymin, dx, dy
            box_height = 1.0
            box_width = 1.0
            population_color_map = self.get_population_colormap(population_data)
            for allele_index, allele in enumerate(sorted_aligned_labeled_repeats):
                sample_id = sorted_sample_ids[allele_index]
                if sample_id.find('-') != -1:
                    sample_id = sample_id[:sample_id.index('-')]
                box_position = [allele_index, 0]
                fcolor = population_color_map.get(sample_id, '#dbdbdb')  # grey default
                ax_bottom.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                                  linewidth=box_line_width + 0.1,
                                                  facecolor=fcolor,
                                                  edgecolor="white"))

            ax_bottom.set_xlim(right=len(aligned_labeled_repeats))
            ax_bottom.set_xticks([x + 0.5 for x in range(len(aligned_labeled_repeats))])
            ax_bottom.set_xticklabels(sorted_sample_ids, ha='center', rotation=xtick_degrees)
            ax_bottom.set_yticks([])
            ax_bottom.set_frame_on(False)
            # ax_bottom.tick_params(axis='x', which='both', length=0)  # if don't want the ticks

            ax_main.set_xlim(right=len(sorted_aligned_labeled_repeats))
            ax_main.set_xticks([])

        if allele_as_row:
            ax_main.set_ylim(top=len(sorted_aligned_labeled_repeats))
            ax_main.set_yticks([x + 0.5 for x in range(len(aligned_labeled_repeats))])
            ax_main.set_yticklabels(sorted_sample_ids[::-1], ha='right', rotation=ylabel_rotation)

            labels = [x for x in range(1, max_repeat_count + 1)]
            ax_main.set_xticks(labels, labels, rotation=xlabel_rotation)
        else:
            ax_main.set_xlim(right=len(sorted_aligned_labeled_repeats))
            ax_main.set_xticks([x + 0.5 for x in range(len(aligned_labeled_repeats))])
            ax_main.set_xticklabels(sorted_sample_ids, ha='center', rotation=xlabel_rotation)

            labels = [y for y in range(1, max_repeat_count + 1)]
            ax_main.set_yticks(labels, labels=labels, rotation=ylabel_rotation)

        ax_main.tick_params(axis='y', which='major', labelsize=ylabel_size)
        ax_main.tick_params(axis='y', which='minor', labelsize=ylabel_size)
        ax_main.tick_params(axis='x', which='major', labelsize=xlabel_size)
        ax_main.tick_params(axis='x', which='minor', labelsize=xlabel_size)

        if hide_xticks:
            ax_main.set_xticks([])
        if hide_yticks:
            ax_main.set_yticks([])

        # Frame
        if frame_on is None:  # Set default
            frame_on = {'top': False, 'bottom': True, 'right': False, 'left': True}

        ax_main.spines['top'].set_visible(frame_on['top'])
        ax_main.spines['right'].set_visible(frame_on['right'])
        ax_main.spines['bottom'].set_visible(frame_on['bottom'])
        ax_main.spines['left'].set_visible(frame_on['left'])

        if output_name is not None:
            if '.' not in output_name:
                fig.savefig(f"{output_name}.pdf", dpi=dpi, bbox_inches='tight')
            else:
                fig.savefig(f"{output_name}", dpi=dpi, bbox_inches='tight')
        else:
            fig.savefig("test_trplot.png", dpi=dpi, bbox_inches='tight')

        plt.close(fig)

    def add_dendrogram(self, fig, aligned_labeled_repeats, sample_ids, symbol_to_motif, hide_clustering):
        """
        Add dendrogram to the figure.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure to add dendrogram to.
        aligned_labeled_repeats : list of list of str
            List of aligned labeled repeats.
        sample_ids : list of str
            List of sample ids.
        symbol_to_motif : dict of str to str
            Dictionary mapping symbols to motifs.
        hide_clustering : bool
            Whether to hide clustering.

        Returns
        -------
        sorted_sample_ids : list of str
            List of sorted sample ids.
        sorted_aligned_labeled_repeats : list of list of str
            List of sorted aligned labeled repeats.
        """
        import scipy.cluster.hierarchy as sch
        from scipy.spatial.distance import squareform
        from trviz.utils import get_distance_matrix
        from trviz.utils import calculate_cost_with_dist_matrix

        distance_matrix = get_distance_matrix(symbol_to_motif)
        dist_mat = np.zeros([len(aligned_labeled_repeats), len(aligned_labeled_repeats)])
        for i in range(len(aligned_labeled_repeats)):
            for j in range(len(aligned_labeled_repeats)):
                tr_1 = aligned_labeled_repeats[i]
                tr_2 = aligned_labeled_repeats[j]
                dist_mat[i, j] = calculate_cost_with_dist_matrix(tr_1, tr_2, distance_matrix, allow_copy_change=True)
        condensed_dist_mat = squareform(dist_mat)

        Y = sch.linkage(condensed_dist_mat, method='single', optimal_ordering=True)
        if hide_clustering:
            Z = sch.dendrogram(Y, orientation='top', get_leaves=True, no_plot=True)
        else:
            ax_top = fig.add_axes([0.125, 0.883, 0.775, 0.2])
            Z = sch.dendrogram(Y, orientation='top', get_leaves=True)
            ax_top.set_xticks([])
            ax_top.set_yticks([])
            ax_top.set_frame_on(False)

        idx2 = Z['leaves']
        sorted_aligned_labeled_repeats = [aligned_labeled_repeats[i] for i in idx2]
        sorted_sample_ids = [sample_ids[i] for i in idx2]

        return sorted_sample_ids, sorted_aligned_labeled_repeats

    @staticmethod
    def get_symbol_to_color_map(alpha, unique_symbol_count, unique_symbols):

        # Color-blind friendly colors (Okabe-Ito)
        if unique_symbol_count <= 7:
            cmap = ListedColormap([(0.902, 0.624, 0),
                                   (0.337, 0.706, 0.914),
                                   (0, 0.620, 0.451),
                                   (0.941, 0.894, 0.259),
                                   (0, 0.447, 0.698),
                                   (0.835, 0.369, 0),
                                   (0.8, 0.475, 0.655)])
        else:
            import distinctipy
            cmap = distinctipy.get_colors(unique_symbol_count, pastel_factor=0.9, rng=777)
            cmap = ListedColormap(cmap)

        symbol_to_color = {r: cmap(i) for i, r in enumerate(list(unique_symbols))}
        return symbol_to_color

    def set_symbol_to_color_map(self, symbol_to_color):
        self.symbol_to_color = symbol_to_color

    def get_population_colormap(self, population_data):
        """ Allowed super populations: {AMR, AFR, EAS, SAS, EUR} """

        # Course et al. 2020 color map
        # popultation_colormap = {'SAS': (117, 213, 253, 100),  # light blue
        #                         'EAS': (219, 219, 219, 100),  # light grey
        #                         'EUR': (211, 250, 121, 100),  # light green+yellow
        #                         'AMR': (252, 211, 119, 100),  # light orange
        #                         'AFR': (255, 137, 215, 100)}  # light pink

        population_colormap = {'SAS': '#fed23f',  # (254, 210, 63)
                               'EAS': '#b5d33d',  # (181, 211, 61)
                               'EUR': '#6ca2ea',  # (108, 162, 234),
                               'AMR': '#442288',  # (68, 34, 136),
                               'AFR': '#eb7d5b',  # (235, 125, 91)
                               }

        population_colormap = {'SAS': '#1982c4',  # blue
                               'EAS': '#8ac926',  # green
                               'EUR': '#ffca3a',  # yellow
                               'AMR': '#6a4c93',  # purple
                               'AFR': '#ff595e',  # red
                               }

        superpopulation_list = ('AMR', 'AFR', 'EAS', 'SAS', 'EUR')

        sample_to_colormap = {}
        with open(population_data, "r") as f:
            for line in f:
                split = line.strip().split("\t")
                sample, population = split[0], split[5]
                sample_to_colormap[sample] = population_colormap.get(population, '#dbdbdb')  # grey color

        return sample_to_colormap


