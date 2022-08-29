import subprocess
import os
import time

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO


class MotifAligner:

    def align(self, sample_ids, labeled_vntrs, vid=None, tool="mafft"):
        motif_aligner = self._get_motif_aligner(tool)
        return motif_aligner(sample_ids, labeled_vntrs, vid)

    def _get_motif_aligner(self, tool):
        if tool == 'mafft':
            return self._align_motifs_with_mafft
        elif tool == 'muscle':
            return self._align_motifs_with_muscle
        elif tool == 'clustalo':
            return self._align_motifs_with_clustalo
        else:
            ValueError(tool)

    @staticmethod
    def _align_motifs_with_muscle(sample_ids, labeled_vntrs, vid):
        muscle_cline = MuscleCommandline('muscle', clwstrict=True)
        data = '\n'.join(['>%s\n' % str(sample_ids[i]) + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        stdout, stderr = muscle_cline(stdin=data)
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_vntrs = [str(aligned.seq) for aligned in alignment]

        return sample_ids, aligned_vntrs  # TODO sample_ids are not correctly sorted

    @staticmethod
    def _align_motifs_with_clustalo(sample_ids, labeled_vntrs, vid):
        clustalo_cline = ClustalOmegaCommandline('clustalo', infile="data.fasta", outfile="test.out",
                                                 force=True,
                                                 clusteringout="cluster.out")

        # Use dist-in (and use pre computed distance)
        # See the clusters - and plot only for those clusters.

        data = '\n'.join(['>%s\n' % str(sample_ids[i]) + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        stdout, stderr = clustalo_cline()
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_vntrs = [str(aligned.seq) for aligned in alignment]

        return sample_ids, aligned_vntrs  # TODO sample_ids are not correctly sorted

    def _align_motifs_with_mafft(self, sample_ids, labeled_vntrs, vid, preserve_order=False):
        temp_input_name = "alignment/alignment_input.fa"
        temp_output_name = "alignment/alignment_output.fa"
        if vid is not None:
            temp_input_name = "alignment/alignment_input_{}.fa".format(vid)
            temp_output_name = "alignment/alignment_output_{}.fa".format(vid)

        data = '\n'.join(['>%s\n' % sample_ids[i] + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        with open(temp_input_name, "w") as f:
            f.write(data)

        # TODO call mafft using pysam wrapper (
        # mafft_exe = "/usr/bin/mafft"
        # mafft_cline = MafftCommandline(mafft_exe, input=temp_input_name)

        if preserve_order:
            os.system("mafft --quiet --text --auto {} > {}".format(temp_input_name, temp_output_name))
        else:
            os.system("mafft --quiet --text --auto --reorder {} > {}".format(temp_input_name, temp_output_name))

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

        if len(aligned_vntrs) == 0:
            raise ValueError(f"No alinged VNTRs in {temp_output_name} file")

        print("sample size", len(sample_ids), len(aligned_vntrs))
        return sample_ids, aligned_vntrs
