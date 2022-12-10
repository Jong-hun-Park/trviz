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

Full documentation is available at [readthedocs](https://trviz.readthedocs.io/)

# Overview of TRviz
<p float="left">
<img src="https://github.com/Jong-hun-Park/trviz/blob/main/examples/figures/TRviz_main_figure.png" width="100%" height="250">
</p>

# Getting Started

## Prerequisite
TRviz requires [MAFFT](https://mafft.cbrc.jp/alignment/software/).

Install the library with pip or from source.
## with pip
```
pip install trviz
```

## from source
```
git clone https://github.com/Jong-hun-Park/trviz.git
cd trviz/
pip install .
```

# Motivation
There have been many approaches to visualize the variations in tandem repeats. 
However, there is no tool available for that.
TRViz automatically decompose tandem repeat sequence into motifs, and align the
decomposed motifs, and finally generate a plot to show the aligned motifs.

## Input
1. Tandem repeat sequences in FASTA format
2. A set of motifs for decomposition

## Output
1. Motif map, a set of motifs detected in the samples and their labels and frequencies
2. Aligned and labeled motifs
3. Plot showing the motif composition of the input sequences
4. Plot mapping color to the motif sequences

# Code samples and examples
TRviz has four modules:
1. Decomposition
2. Encoding
3. Alignment
4. Visualization

See full documentation at [readthedocs]()

## Generating a TR plot

```python
from trviz.main import TandemRepeatVizWorker
from trviz.utils import get_sample_and_sequence_from_fasta

tr_visualizer = TandemRepeatVizWorker()
sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(fasta_file_path)
tr_id = "CACNA1C"
motifs = ['GACCCTGACCTGACTAGTTTACAATCACAC']

tr_visualizer.generate_trplot(tr_id, sample_ids, tr_sequences, motifs)
``` 

## Motif Decomposition
```python
from trviz.decomposer import Decomposer

tr_decomposer = Decomposer()
tr_sequence = "ACCTTGACCTTGACCTTGACCTTG"
motifs = ["ACCTTG"]
tr_decomposer.decompose(tr_sequence, motifs)
# >>> ["ACCTTG", "ACCTTG", "ACCTTG", "ACCTTG"]
``` 

# Contact Us
Please submit an issue on the [TRviz github](https://github.com/Jong-hun-Park/trviz/issues)
