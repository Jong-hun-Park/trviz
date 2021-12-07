try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


def align_motifs(labeled_vntrs):
    muscle_cline = MuscleCommandline('muscle', clwstrict=True)
    data = '\n'.join(['>%s\n' % str(i) + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
    stdout, stderr = muscle_cline(stdin=data)
    alignment = AlignIO.read(StringIO(stdout), "clustal")
    aligned_motifs = [str(aligned.seq) for aligned in alignment]

    return aligned_motifs
