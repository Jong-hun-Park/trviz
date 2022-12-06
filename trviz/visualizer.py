from typing import List, Tuple, Dict
import re
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

        fig, ax = plt.subplots()
        if figure_size is not None:
            fig, ax = plt.subplots(figsize=figure_size)
        ax.set_aspect('equal')

        max_motif_length = len(max(symbol_to_motif.values(), key=len))
        # TODO: sort by the frequency
        # TODO: figure and font size
        yticklabels = []
        y_base_position = 0
        for (symbol, color) in symbol_to_color.items():
            if symbol == '-':
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
        ax.set_yticks([y + 0.5 for y in range(len(symbol_to_color) - 1)])  # minus 1 for gap symbol, '-'
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
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)

        fig.tight_layout()
        fig.savefig(file_name, dpi=dpi, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

    def trplot(self,
               aligned_labeled_repeats: List[str],
               sample_ids: List[str] = None,
               figure_size: Tuple[int, int] = None,
               output_name: str = None,
               dpi: int = 500,
               alpha: float = 0.6,
               box_line_width: float = 0,
               xtick_degrees: int = 90,
               hide_xticks: bool = False,
               hide_yticks: bool = False,
               symbol_to_motif: Dict = None,
               sort_by_clustering: bool = True,
               hide_dendrogram: bool = False,
               population_data = None,
               sample_as_row = False,
               private_motif_color: str = 'black',
               frame_on = {'top': False, 'bottom': True, 'right': False, 'left': True},
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
        :param xtick_degrees: xtick degree (default is 90)
        :param hide_xticks: if true, hide xticks
        :param hide_yticks: if true, hide yticks
        :param private_motif_color: the color for private motifs. Default is black
        :param debug: if true, print verbse information.
        """

        fig, ax_main = plt.subplots(figsize=figure_size)  # width and height in inch

        max_repeat_count = len(aligned_labeled_repeats[0])
        if figure_size is None:
            w = len(sample_ids) // 10 + 5 if len(sample_ids) > 50 else 5
            h = max_repeat_count // 10 + max_repeat_count if max_repeat_count > 50 else max_repeat_count // 5 + 2
            if h * dpi > 2**16:
                h = 2**15 // dpi * 0.75  # "Weight and Height must be less than 2^16"
            fig, ax_main = plt.subplots(figsize=(w, h))

        # Add clustering
        if sort_by_clustering:
            if debug:
                print("Sort by clustering")
            sorted_sample_ids, sorted_aligned_labeled_repeats = self.add_dendrogram(fig,
                                                                                    aligned_labeled_repeats,
                                                                                    sample_ids,
                                                                                    symbol_to_motif,
                                                                                    hide_dendrogram)
        else:
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
            box_position = [allele_index, 0]
            for box_index, symbol in enumerate(allele):
                box_position[1] = box_height * box_index
                fcolor = symbol_to_color[symbol]
                if symbol == '-':  # For gaps, color them as white blocks
                    fcolor = (1, 1, 1, 1)
                if symbol == PRIVATE_MOTIF_LABEL:
                    fcolor = private_motif_color
                ax_main.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                           linewidth=box_line_width + 0.1,
                                           facecolor=fcolor,
                                           edgecolor="white"))

        if population_data is not None:
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
        else:
            ax_main.set_xlim(right=len(sorted_aligned_labeled_repeats))
            ax_main.set_xticks([x + 0.5 for x in range(len(aligned_labeled_repeats))])
            ax_main.set_xticklabels(sorted_sample_ids, ha='center', rotation=xtick_degrees)

        ax_main.set_yticks([y for y in range(1, max_repeat_count + 1)])
        sample_as_row = True
        if sample_as_row:
            # y-tick right and rotate
            ax_main.set_yticklabels(([y for y in range(1, max_repeat_count + 1)]), rotation=90, ha='center')
            ax_main.yaxis.tick_right()
            frame_on['right'] = True
            frame_on['left'] = False
        else:
            ax_main.set_yticklabels(([y for y in range(1, max_repeat_count + 1)]))

        ax_main.tick_params(axis='y', which='major', pad=6.5)
        ax_main.tick_params(axis='both', which='major', labelsize=8)
        ax_main.tick_params(axis='y', which='major', labelsize=9)
        ax_main.tick_params(axis='y', which='minor', labelsize=9)
        ax_main.tick_params(axis='x', which='major', labelsize=9)
        ax_main.tick_params(axis='x', which='minor', labelsize=9)

        if hide_xticks:
            ax_main.set_xticks([])
        if hide_yticks:
            ax_main.set_yticks([])

        # Frame
        ax_main.spines['top'].set_visible(frame_on['top'])
        ax_main.spines['right'].set_visible(frame_on['right'])
        ax_main.spines['bottom'].set_visible(frame_on['bottom'])
        ax_main.spines['left'].set_visible(frame_on['left'])

        if output_name is not None:
            # fig.set_size_inches(w, h)  # Because the figure size is changed after adding axes
            # if re.match(".*", tweet):
            #     fig.savefig(f"{output_name}", dpi=dpi, bbox_inches='tight')
            # else:
            fig.savefig(f"{output_name}.pdf", dpi=dpi, bbox_inches='tight')
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
            cmap = plt.cm.get_cmap('hsv', unique_symbol_count + 1)
            # Setting alpha values
            temp_cmap = cmap(np.arange(cmap.N))
            temp_cmap[:, -1] = alpha
            # Update color map
            cmap = ListedColormap(temp_cmap)

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


