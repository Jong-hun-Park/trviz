
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap


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


def get_unique_repeats(aligned_repeats):
    unique_repeats = set()
    for rs in aligned_repeats:
        for r in rs:
            unique_repeats.add(r)

    return sorted(list(unique_repeats))


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


def trplot(aligned_repeats,
           figure_size=None,
           outfolder=None, outname=None,
           dpi=700,
           alpha=0.5,
           xticks=None,
           xtick_degrees=90,
           hide_yticks=False,
           debug=False):
    """

    :type data: object
    """
    # Load the data as dataframe
    # input is
    # ['BBB----', 'BBBAB--', 'BBCABBB']
    # Convert alphabet to number

    fig = plt.figure()
    if figure_size is not None:
        fig = plt.figure(figsize=figure_size)  # width and height in inch

    max_repeat = len(aligned_repeats[0])
    if debug:
        print("Max repeats: {}".format(max_repeat))

    ax = plt.subplot(1, 1, 1, aspect=1.5, autoscale_on=False, frameon=False,
                     xticks=[x for x in range(len(aligned_repeats) + 1)],
                     yticks=[y for y in range(max_repeat + 1)])

    unique_repeats = get_unique_repeats(aligned_repeats)
    number_of_colors = len(unique_repeats)

    # cmap = plt.cm.get_cmap("tab20")  #XXX Note that the number of colors maybe not enough!
    cmap = plt.cm.get_cmap('hsv', number_of_colors)

    # Setting alpha values
    temp_cmap = cmap(np.arange(cmap.N))
    temp_cmap[:, -1] = alpha
    # New colormap
    cmap = ListedColormap(temp_cmap)

    # colors = cmap(len(unique_repeats))
    motif_color_map = {r: cmap(i) for i, r in enumerate(list(unique_repeats))}

    allele_count = len(aligned_repeats)
    box_height = max_repeat/len(aligned_repeats[0])
    box_width = 1.0
    # print("box height", box_height)
    # print("box width", box_width)

    for allele_index, allele in enumerate(aligned_repeats):
        if debug:
            print("Allele index", allele_index)
        box_position = [allele_index, 0]
        for box_index, repeat in enumerate(allele):
            box_position[1] = box_height * box_index
            if debug:
                print("Drawing box in position", box_position)
            fcolor = motif_color_map[repeat]
            if repeat == '-':  # For gaps, color them as white blocks
                fcolor = (1, 1, 1, 1)
            ax.add_patch(plt.Rectangle(box_position, box_width, box_height,
                                       linewidth=0.1,
                                       facecolor=fcolor,
                                       edgecolor="white"))

    if xticks is not None:
        plt.xticks([x+0.5 for x in range(len(xticks))], xticks,
                   fontsize=8,
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
    if outname is not None and outfolder is not None:
        plt.savefig("{}/{}.png".format(outfolder, outname), dpi=dpi)
    else:
        plt.savefig("test_vntr_plot.png", dpi= dpi)

    plot_motif_color_map(motif_color_map, file_name="{}/motif_color_map_{}.png".format(outfolder, outname))
