Code samples
============

## Decomposition
```python
from trviz.decomposer import Decomposer

tr_decomposer = Decomposer()
tr_sequence = "ACCTTGACCTTGACCTTGACCTTG"
motifs = ["ACCTTG"]
tr_decomposer.decompose(tr_sequence, motifs)
# >>> ['ACCTTG', 'ACCTTG', 'ACCTTG', 'ACCTTG']
```

## Encoding
```python
from trviz.motif_encoder import MotifEncoder

motif_encoder = MotifEncoder()
decomposed_vntrs = [['ACCTTG', 'ACCTTG', 'ACCTTC'],
                    ['ACCTTG', 'ACCTTG', 'ACCTTG', 'ACCTTC'],
                    ['ACCTTG', 'ACCTTG', 'ACCTTG', 'ACCTTC', 'ACCTTC'],
                   ]
motif_encoder.encode(decomposed_vntrs, motif_map_file="VNTR_motif_map.txt")
# >>> ['aab', 'aaab', 'aaabb']
```

## Alignment
```python
from trviz.motif_aligner import MotifAligner

motif_aligner = MotifAligner()
motif_aligner.align(sample_ids = ['sample1', 'sample2', 'sample3'],
                    encoded_vntrs = ['aab', 'aaab', 'aaabb'],
                    vid = 'test')
# >>> (['sample1', 'sample2', 'sample3'], ['-aab-', 'aaab-', 'aaabb'])
```

## Visualization
```python
from trviz.visualizer import TandemRepeatVisualizer

sample_ids = ['sample1', 'sample2', 'sample3']
aligned_vntrs = ['-aab-', 'aaab-', 'aaabb']
visualizer = TandemRepeatVisualizer()
# This will generate a TR plot as a file, "test_tr_plot.png"
visualizer.trplot(aligned_vntrs, output_name="test_tr_plot.png", sample_ids=sample_ids)

```
