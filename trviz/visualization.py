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

        # fig = plt.figure(figsize=(8, 5))  # width and height in inch
        ax = plt.subplot(1, 1, 1, aspect=1, autoscale_on=True, frameon=False,
                         xticks=[y for y in range(len(motif_color_map) + 1)],
                         yticks=[y for y in range(len(motif_color_map) + 1)])

        for i, (motif, color) in enumerate(motif_color_map.items()):
            print(motif)
            print(color)
            print(i)
            box_position = [0, i]
            box_width = 0.5
            box_height = 0.5
            ax.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                       linewidth=0,
                                       facecolor=color,
                                       edgecolor="white"))
            text_position = [box_position[0] + 2.5 * box_width / 2, box_position[1] + box_height / 2]
            ax.text(x=text_position[0], y=text_position[1], s=motif, fontsize=20)

        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        plt.savefig(file_name)

    def plot(self,
             aligned_labeled_repeats,
             figure_size=None,
             output_name=None,
             dpi=500,
             alpha=0.5,
             sample_ids=None,
             xtick_degrees=90,
             hide_yticks=False,
             private_motif_color='black',
             debug=False):

        """ Generate tandem repeat plot """

        fig = plt.figure()
        if figure_size is not None:
            fig = plt.figure(figsize=figure_size)  # width and height in inch

        max_repeat = len(aligned_labeled_repeats[0])
        x_y_ratio = max_repeat / len(sample_ids)
        if debug:
            print("Max repeats: {}".format(max_repeat))
            print("x vs y ratio", x_y_ratio)

        ax = fig.add_subplot(1, 1, 1,
                             # aspect=1.5,
                             aspect=1/(x_y_ratio)*4,
                             autoscale_on=False,
                             frameon=False,
                             xticks=[x for x in range(len(aligned_labeled_repeats) + 1)],
                             yticks=[y for y in range(max_repeat + 1)])

        unique_labels = self._get_unique_labels(aligned_labeled_repeats)
        unique_label_count = len(unique_labels)

        if unique_label_count < 20:
            cmap = plt.cm.get_cmap("tab20")
        else:
            cmap = plt.cm.get_cmap('hsv', unique_label_count)

        # Setting alpha values
        temp_cmap = cmap(np.arange(cmap.N))
        temp_cmap[:, -1] = alpha
        # Update color map
        cmap = ListedColormap(temp_cmap)

        # colors = cmap(len(unique_labels))
        label_to_color = {r: cmap(i) for i, r in enumerate(list(unique_labels))}

        box_height = max_repeat/len(aligned_labeled_repeats[0])
        box_width = 1.0

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
                                           linewidth=0.1,
                                           facecolor=fcolor,
                                           edgecolor="white"))

        if sample_ids is not None:
            plt.xticks([x + 0.5 for x in range(len(sample_ids))],
                       sample_ids,
                       fontsize=3,
                       rotation=xtick_degrees)
        else:
            plt.xticks([])

        if hide_yticks:
            plt.yticks([])

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=True)  # labels along the bottom edge are off

        plt.tight_layout()

        if output_name is not None:
            plt.savefig("{}.png".format(output_name), dpi=dpi)
        else:
            plt.savefig("test_vntr_plot.png", dpi=dpi)
