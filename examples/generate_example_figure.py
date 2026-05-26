"""Generate a tandem repeat plot for VNTR in SORL1 gene.

Demonstrates two patterns:

1. The original "save to file" usage (uncomment to run).
2. The v1.4.0 fig/ax return for matplotlib-idiomatic post-styling.
"""

from trviz.main import TandemRepeatVizWorker
from trviz.utils import get_sample_and_sequence_from_fasta

fasta_file = "data/SORL1_tr_sequence.fa"
sample_ids, tandem_repeat_sequences = get_sample_and_sequence_from_fasta(fasta_file)
tandem_repeat_id = "SORL1"
motifs = ["GCTTCATCTCCTCCTCCTCACCTCCTGCTGTGGTGCACAGATACCTATAGGCAG"]

tr_visualizer = TandemRepeatVizWorker()

# Pattern 1: write the plot directly to disk (old behavior, still the default).
tr_visualizer.generate_trplot(
    tandem_repeat_id,
    sample_ids,
    tandem_repeat_sequences,
    motifs,
    figure_size=(10, 13),
)

# Pattern 2 (v1.4.0+): get back a (fig, ax) tuple and apply your own matplotlib styling.
# Uncomment to try:
#
# fig, ax = tr_visualizer.generate_trplot(
#     tandem_repeat_id,
#     sample_ids,
#     tandem_repeat_sequences,
#     motifs,
#     figure_size=(10, 13),
#     save=False,           # skip the built-in savefig
# )
# ax.set_title("SORL1 VNTR composition", fontsize=14)
# ax.tick_params(axis="x", labelsize=12)
# fig.savefig("SORL1_custom.pdf", dpi=300, bbox_inches="tight")
