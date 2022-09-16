from trviz.main import TandemRepeatVizWorker
from trviz.utils import read_fasta


""" Generate a Tandem Repeat plot for VNTR in SORL1 gene """

fasta_file = "data/SORL1_tr_sequence.fa"
sample_ids, tandem_repeat_sequences = read_fasta(fasta_file)
tandem_repeat_id = "SORL1"
motifs = ['GCTTCATCTCCTCCTCCTCACCTCCTGCTGTGGTGCACAGATACCTATAGGCAG']

tr_visualizer = TandemRepeatVizWorker()
tr_visualizer.generate_tr_plot(tandem_repeat_id,
                               sample_ids,
                               tandem_repeat_sequences,
                               motifs,
                               figure_size=(10, 13),
                               )
