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
## Prerequisite
> **Note**
> Before getting started, ensure you have [MAFFT](https://mafft.cbrc.jp/alignment/software/).
> The current version is tested with MAFFT v7.505.

## Step 1: Install TRviz

```bash
# Install with pip
pip install trviz
```
or
```bash
# Install from source
git clone https://github.com/Jong-hun-Park/trviz.git
cd trviz/
pip install .
```

## Step 2: Run Your First Analysis
Check out our [Jupyter Notebook](https://github.com/Jong-hun-Park/trviz/blob/main/examples/sample_code.ipynb) for code examples and start visualizing tandem repeat sequence right away!

# Features Overview
![](https://github.com/Jong-hun-Park/trviz/blob/main/examples/figures/TRviz_main_figure.png?raw=true)


# Installation
Use pip for a quick installation, or install from source for more control.
```bash
pip install trviz
````
or 
```bash
# Install from source
git clone https://github.com/Jong-hun-Park/trviz.git
cd trviz/
pip install .
```

# Input and Output
## Input
1. Tandem repeat sequences (alleles)
2. A set of motifs for decomposition

## Output
1. A plot showing the motif composition of the input sequences (pdf by default)
2. A plot mapping color to motif (pdf by default)
3. Aligned and labeled motifs (text file)
4. Motif map, a set of motifs detected in the samples and their labels and frequencies (text file) 

For more detailed descriptions, please see full documentation at [readthedocs](https://trviz.readthedocs.io/en/latest/)

# Code examples
## Generating a plot

```python
from trviz.main import TandemRepeatVizWorker
from trviz.utils import get_sample_and_sequence_from_fasta

tr_visualizer = TandemRepeatVizWorker()
sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(fasta_file_path)
tr_id = "CACNA1C"
motifs = ['GACCCTGACCTGACTAGTTTACAATCACAC']

tr_visualizer.generate_trplot(tr_id, sample_ids, tr_sequences, motifs)
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

# Citation:
Jonghun Park, Eli Kaufman, Paul N Valdmanis, Vineet Bafna, 
[TRviz: a Python library for decomposing and visualizing tandem repeat sequences](https://doi.org/10.1093/bioadv/vbad058), Bioinformatics Advances, Volume 3, Issue 1, 2023, vbad058


# Contribute
Your feedback is valuable! If you encounter any issues during installation or usage, please submit them in the [TRviz GitHub Issues](https://github.com/Jong-hun-Park/trviz/issues).
