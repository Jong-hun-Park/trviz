try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO


def align_motifs(labeled_vntrs, method="muscle"):
    if method == "muscle":
        muscle_cline = MuscleCommandline('muscle', clwstrict=True)
        data = '\n'.join(['>%s\n' % str(i) + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        stdout, stderr = muscle_cline(stdin=data)
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_motifs = [str(aligned.seq) for aligned in alignment]

    elif method == "clustalo":
        clustalo_cline = ClustalOmegaCommandline('clustalo', infile="data.fasta", outfile="test.out",
                                                 force=True,
                                                 clusteringout="cluster.out")

        # Use dist-in (and use pre computed distance)
        # See the clusters - and plot only for those clusters.

        data = '\n'.join(['>%s\n' % str(i) + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        stdout, stderr = clustalo_cline()
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_motifs = [str(aligned.seq) for aligned in alignment]

    elif method == "mafft":

        data = '\n'.join(['>%s\n' % str(i) + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        with open("test_input.fasta", "w") as f:
            f.write(data)
        mafft_cline = MafftCommandline(input="test_input.fasta")
        print(mafft_cline)
        stdout, stderr = mafft_cline()
        print(stdout)
        import subprocess
        # import os
        # stdout = os.system("mafft --anysymbol < {}".format(data))
        # print(stdout)
        # stdout = subprocess.run(["mafft", "--anysymbol", "<", data], stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)

        alignment = AlignIO.read(StringIO(stdout), "fasta")
        aligned_motifs = [str(aligned.seq) for aligned in alignment]

    return aligned_motifs
