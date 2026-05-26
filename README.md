```
  _______ _______      _______ ______
 |__   __|  __ \ \    / /_   _|___  /
    | |  | |__) \ \  / /  | |    / /
    | |  |  _  / \ \/ /   | |   / /
    | |  | | \ \  \  /   _| |_ / /__
    |_|  |_|  \_\  \/   |_____/_____|

```
TRviz is a python library for analyzing tandem repeat sequences. TRviz includes modules for
decomposing, encoding, aligning, and visualizing tandem repeat sequences.

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://trviz.readthedocs.io/)

# Quick Start

## Prerequisites

- **Python 3.10+** (3.10, 3.11, or 3.12 are tested)
- **MAFFT** (for multiple sequence alignment) — install from [MAFFT website](https://mafft.cbrc.jp/alignment/software/). Tested with v7.505.

## Install

```bash
pip install trviz
```

Or from source:

```bash
git clone https://github.com/Jong-hun-Park/trviz.git
cd trviz/
pip install .
```

## Run your first analysis

Check out our [Jupyter notebook](https://github.com/Jong-hun-Park/trviz/blob/main/examples/sample_code.ipynb) for code examples.

# Features Overview

![](https://github.com/Jong-hun-Park/trviz/blob/main/examples/figures/TRviz_main_figure.png?raw=true)

# Input and output

## Input
1. Tandem repeat sequences (alleles)
2. A set of motifs for decomposition

## Output
1. A plot showing the motif composition of the input sequences (PDF by default)
2. A plot mapping color to motif (PDF by default)
3. Aligned and labeled motifs (text file)
4. Motif map: a set of motifs detected in the samples and their labels and frequencies (text file)

For more detailed descriptions, see the full documentation at [readthedocs](https://trviz.readthedocs.io/en/latest/).

# Code examples

## Generate a plot (saves to disk)

```python
from trviz.main import TandemRepeatVizWorker
from trviz.utils import get_sample_and_sequence_from_fasta

tr_visualizer = TandemRepeatVizWorker()
sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(fasta_file_path)
tr_id = "CACNA1C"
motifs = ['GACCCTGACCTGACTAGTTTACAATCACAC']

tr_visualizer.generate_trplot(tr_id, sample_ids, tr_sequences, motifs)
```

## Customize the plot via matplotlib (v1.4.0+)

`generate_trplot()` and the underlying `trplot()` / `plot_motif_color_map()` now return `(fig, ax)`,
so you can apply any matplotlib styling — fonts, ticks, titles, layout — before saving the figure
yourself.

```python
fig, ax = tr_visualizer.generate_trplot(
    tr_id, sample_ids, tr_sequences, motifs,
    save=False,                       # skip the built-in savefig
)

# Apply any matplotlib styling you want:
ax.set_title("CACNA1C VNTR composition", fontsize=14)
ax.tick_params(axis="x", labelsize=12)
ax.set_xlabel("Repeat unit", fontsize=12)

fig.savefig("CACNA1C.pdf", dpi=300, bbox_inches="tight")
```

You can also draw the plot into your own `ax` (useful for embedding in a multi-panel figure):

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
tr_visualizer.visualizer.trplot(
    aligned_labeled_repeats=aligned_trs,
    sample_ids=sorted_sample_ids,
    symbol_to_motif=symbol_to_motif,
    ax=axes[0],                       # draw into your own axes
    sort_by_clustering=False,
)
# ... draw something else into axes[1] ...
fig.savefig("multi_panel.pdf")
```

## Motif decomposition

```python
from trviz.decomposer import Decomposer

tr_decomposer = Decomposer()
tr_sequence = "ACCTTGACCTTGACCTTGACCTTG"
motifs = ["ACCTTG"]
tr_decomposer.decompose(tr_sequence, motifs)
# >>> ["ACCTTG", "ACCTTG", "ACCTTG", "ACCTTG"]
```

# Citation

Jonghun Park, Eli Kaufman, Paul N Valdmanis, Vineet Bafna,
[TRviz: a Python library for decomposing and visualizing tandem repeat sequences](https://doi.org/10.1093/bioadv/vbad058),
Bioinformatics Advances, Volume 3, Issue 1, 2023, vbad058.

# Contribute

Your feedback is valuable! If you encounter any issues during installation or usage, please submit
them in the [TRviz GitHub Issues](https://github.com/Jong-hun-Park/trviz/issues).
