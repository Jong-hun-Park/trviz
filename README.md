# TRviz
A tool to decompose, align, and visualize tandem repeat sequences.

## Examples
<p float="left">
<img src="https://github.com/Jong-hun-Park/TandemRepeatViz/blob/main/examples/example_figure1_VPS53.png" width="32%" height="300">
<img src="https://github.com/Jong-hun-Park/TandemRepeatViz/blob/main/examples/example_figure2_SORL1.png" width="32%" height="300">
<img src="https://github.com/Jong-hun-Park/TandemRepeatViz/blob/main/examples/example_figure3_CACNA1.png" width="32%" height="300">
</p>

## Installation
TBD

## Motivation
There have been many approaches to visualize the variations in tandem repeats. 
However, there is no tool available for that.
TRViz automatically decompose tandem repeat sequence into motifs, and align the
decomposed motifs, and finally generate a plot to show the aligned motifs.

### Input
1. Tandem repeat sequences in FASTA format
2. (Optional) Motifs for decomposition

### Output
1. Motif map, a set of motifs detected in the samples and their labels and frequencies
2. Aligned labeled motifs
3. Plot


## Usage
TRViz has three modules:
1. Decomposition
2. Alignment
3. Visualization

```python
from trviz.main import TandemRepeatVizWorker
from trviz.utils import read_fasta
tr_visualizer = TandemRepeatVizWorker()
samples, tr_sequences = read_fasta(input_fasta)
vntr_id = "CACNA1C"
motifs = ['GACCCTGACCTGACTAGTTTACAATCACAC']

tr_visualizer.generate_tr_plot(tr_sequences,
                               motifs,
                               vntr_id,
                               samples,
                               )
``` 
