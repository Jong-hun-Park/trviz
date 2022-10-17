from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from trviz.utils import PRIVATE_MOTIF_LABEL


class TandemRepeatVisualizer:

    def __init__(self):
        pass

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
    def plot_motif_color_map(motif_color_map, file_name):
        """ Generate a plot for motif color map """
        # fig = plt.figure(figsize=(8, 5))  # width and height in inch
        ax = plt.subplot(1, 1, 1, aspect=1, autoscale_on=True, frameon=False,
                         xticks=[y for y in range(len(motif_color_map) + 1)],
                         yticks=[y for y in range(len(motif_color_map) + 1)])

        for i, (motif, color) in enumerate(motif_color_map.items()):
            box_position = [0, i]
            box_width = 0.5
            box_height = 0.5
            ax.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                       linewidth=0,
                                       facecolor=color,
                                       edgecolor="white"))
            text_position = [box_position[0] + 2.5 * box_width / 2, box_position[1] + box_height / 2 - 0.1]
            ax.text(x=text_position[0], y=text_position[1], s=motif, fontsize=20)

        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        plt.savefig(file_name)
        plt.close()

    def plot(self,
             aligned_labeled_repeats: List[str],
             figure_size: Tuple[int, int] = None,
             output_name: str = None,
             dpi: int = 500,
             alpha: float = 0.5,
             sample_ids: List[str] = None,
             box_line_width: float = 0,
             xtick_degrees: int = 90,
             hide_xticks: bool = False,
             hide_yticks: bool = False,
             private_motif_color: str = 'black',
             debug: bool = False
             ):
        """
        Generate a plot showing the variations in tandem repeat sequences.
        A distinct color is assigned to each public motif. For private motifs, the same color (black by default) is assigned.

        :param aligned_labeled_repeats: aligned and encoded tandem repeat sequences.
        :param figure_size: figure size
        :param output_name: output file name
        :param dpi: DPI for the plot
        :param alpha: alpha value for the plot
        :param sample_ids: sample IDs
        :param box_line_width: line width for box edges
        :param xtick_degrees: xtick degree (default is 90)
        :param hide_xticks: if true, hide xticks
        :param hide_yticks: if true, hide yticks
        :param private_motif_color: the color for private motifs. Default is black
        :param debug: if true, print verbse information.
        """

        max_repeat = len(aligned_labeled_repeats[0])
        x_y_ratio = max_repeat / len(sample_ids)
        if debug:
            print("Max repeats: {}".format(max_repeat))
            print("x vs y ratio", x_y_ratio)

        unique_labels = self._get_unique_labels(aligned_labeled_repeats)
        unique_label_count = len(unique_labels)

        # if unique_label_count < 20:
        #     cmap = plt.cm.get_cmap("tab20")
        # else:
        #     cmap = plt.cm.get_cmap('hsv', unique_label_count)
        cmap = plt.cm.get_cmap('hsv', unique_label_count)

        # Setting alpha values
        temp_cmap = cmap(np.arange(cmap.N))
        temp_cmap[:, -1] = alpha
        # Update color map
        cmap = ListedColormap(temp_cmap)

        # colors = cmap(len(unique_labels))
        label_to_color = {r: cmap(i) for i, r in enumerate(list(unique_labels))}

        # self.plot_motif_color_map(label_to_color, output_name + "_color_map.png")

        fig, ax = plt.subplots(figsize=figure_size)  # width and height in inch
        if figure_size is None:
            w = len(sample_ids) // 10 + 3 if len(sample_ids) > 50 else 5
            h = max_repeat // 10 + 5 if max_repeat > 50 else max_repeat // 5 + 2
            if h * dpi > 2**16:
                h = 2**15 // dpi  # "Weight and Height must be less than 2^16"
            fig, ax = plt.subplots(figsize=(w, h))

        box_height = max_repeat/len(aligned_labeled_repeats[0])
        box_width = 1.0

        ax.set_xticks([x + box_width/2 for x in range(len(aligned_labeled_repeats))])
        ax.set_yticks([y for y in range(max_repeat + 1)])
        ax.set_ylabel("Motif counts")

        for allele_index, allele in enumerate(aligned_labeled_repeats):
            box_position = [allele_index, 0]
            for box_index, label in enumerate(allele):
                box_position[1] = box_height * box_index
                fcolor = label_to_color[label]
                if label == '-':  # For gaps, color them as white blocks
                    fcolor = (1, 1, 1, 1)
                if label == PRIVATE_MOTIF_LABEL:
                    fcolor = private_motif_color
                ax.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                           linewidth=box_line_width,
                                           facecolor=fcolor,
                                           edgecolor="white"))

        if not hide_xticks:
            ax.set_xticklabels(sample_ids,
                               rotation=xtick_degrees)

        else:
            ax.set_xticklabels.xticks([])

        if hide_yticks:
            ax.set_yticklabels.yticks([])

        ax.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=True)  # labels along the bottom edge are off

        fig.tight_layout()

        if output_name is not None:
            plt.savefig("{}.png".format(output_name), dpi=dpi)
        else:
            plt.savefig("test_vntr_plot.png", dpi=dpi)

        plt.close(fig)
