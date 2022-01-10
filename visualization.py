
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns


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


def trplot(data):
    """

    :type data: object
    """
    # Load the data as dataframe
    # input is
    # ['BBB----', 'BBBAB--', 'BBCABBB']
    # Convert alphabet to number

    encoded_data = encode_tr_sequence(data)
    print(encoded_data)

    max_ru_length = len(encoded_data[0])
    df_dict = {}
    for i in range(max_ru_length):
        df_dict[str(i+1)] = list()
        for d in encoded_data:
            df_dict[str(i+1)].append(d[i])

    print(df_dict)
    df = pd.DataFrame.from_dict(df_dict)
    print(df)
    # print(df[[1, 2, 3]])

    df = pd.read_csv("test_dict.csv")
    df.set_index('sample', inplace=True)
    print("DF")
    print(df)


    fig = plt.figure()

    num_colors = df.to_numpy().max()
    print("num colors", num_colors)
    cmap = ListedColormap([(1.0, 1.0, 1.0, 1.0)] + sns.color_palette("plasma", n_colors=num_colors))

    g = sns.clustermap(df[[str(i+1) for i in range(len(df.columns))]],
                       col_cluster=False,
                       linewidth=0.2,
                       # figsize=(len(df.columns), len(df.index)/10),
                       yticklabels=False,
                       cbar_kws={'ticks': [i for i in range(num_colors+1)]},
                       vmin=0,
                       vmax=num_colors,
                       cmap=cmap)

    g.cax.set_visible(False)  # No colorbar
    ax = g.ax_heatmap

    ax.set_xlabel("Repeat Counts", fontsize=20)
    ax.set_xlabel("Alleles", fontsize=20)

    # plt.show()
    plt.tight_layout()
    plt.savefig("test.png")
    # plt.close()
