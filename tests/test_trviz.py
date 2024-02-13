from trviz.main import TandemRepeatVizWorker

tr_visualizer = TandemRepeatVizWorker()
tr_id = "test"
tr_seqs = ['ACTACTACTACTACT',
           'ACTACTACTACTACTACTACTACTACTACTACTACTACT',
           'ACTACTACTACTACT',
           'ACTACTACTACTACTACTACTACTACTACTACTACT',
           'ACTACTACTACTACTACTACTACTACTACTACT',
           'ACTACTACTACTACTACTACTACTACTACT',
           'ACTACTACTACTACTACTACTACTACT',
           'ACTACTACTACTACTACTACTACT',
           'ACTACTACTACTACTACT',
           'ACTACTACTACTACTACTACT',
           ]
tr_seqs = tr_seqs * 2

sample_ids = [f"sample{x}" for x in range(1, len(tr_seqs) + 1)]
motifs = ['ACT']
tr_visualizer.generate_trplot(tr_id, sample_ids, tr_seqs, motifs, hide_dendrogram=False, skip_alignment=True)
