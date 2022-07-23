try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO


def align_motifs(sample_ids, labeled_vntrs, vid=None, method="muscle"):
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
        temp_input_name = "temp/temp_input.fa"
        temp_output_name = "temp/temp_output.fa"
        if vid is not None:
            temp_input_name = "temp/temp_input_{}.fa".format(vid)
            temp_output_name = "temp/temp_output_{}.fa".format(vid)

        data = '\n'.join(['>%s\n' % sample_ids[i] + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        with open(temp_input_name, "w") as f:
            f.write(data)

        # mafft_cline = MafftCommandline(input="test_input.fasta")
        # print(mafft_cline)
        # stdout, stderr = mafft_cline()
        # print(stdout)

        import subprocess
        import os
        stdout = os.system("mafft --text --auto {} > {}".format(temp_input_name, temp_output_name))

        # import subprocess
        # align_out = open(temp_output_name, "w")
        # p = subprocess.run(['mafft', '--anysymbol', '--auto', temp_input_name],
        #                    shell=True,
        #                    stdout=align_out,
        #                    stderr=subprocess.STDOUT)
        # print(p)

        import time
        print("sleeping...")
        time.sleep(3)
        # stdout = subprocess.run(["mafft", "--anysymbol", "temp_input.fasta"], stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
        # print(stdout)

        # alignment = AlignIO.read(StringIO(stdout), "fasta")
        # aligned_motifs = [str(aligned.seq) for aligned in alignment]

        aligned_vntrs = []
        sample_ids = []
        tr_seq = None
        with open(temp_output_name, "r") as f:
            for line in f:
                if line.startswith(">"):
                    sample_ids.append(line.strip()[1:])
                    if tr_seq is not None:
                        aligned_vntrs.append(tr_seq)
                    tr_seq = ""
                else:
                    tr_seq += line.strip()

        if len(tr_seq) > 0:
            aligned_vntrs.append(tr_seq)

    print("sample size", len(sample_ids), len(aligned_vntrs))
    return sample_ids, aligned_vntrs
