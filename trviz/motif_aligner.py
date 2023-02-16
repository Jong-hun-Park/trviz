import subprocess
import os
from typing import Tuple
from typing import List
from typing import Dict

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

import shutil


class MotifAligner:
    def align(self,
              sample_ids: List[str],
              encoded_vntrs: List[str],
              vid: str = None,
              score_matrix: Dict[Dict, int] = None,
              output_dir: str = "./",
              tool: str = "mafft",
              ) -> Tuple[List, List]:
        """
        Align encoded VNTRs using multiple sequence alignment tools. Default tool is MAFFT.

        :param sample_ids: sample ids
        :param encoded_vntrs: encoded tandem repeats
        :param tool: the tool name for multiple sequence alignment (options: MAFFT (default))
        :param score_matrix: defines the score matrix between all pair of motifs
        :param output_dir: base directory for output file
        :param vid: ID for the tandem repeat
        """
        motif_aligner = self._get_motif_aligner(tool)
        return motif_aligner(sample_ids, encoded_vntrs, vid, score_matrix, output_dir)

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
    def _align_motifs_with_muscle(sample_ids, labeled_vntrs, vid, score_matrix, output_dir):
        muscle_cline = MuscleCommandline('muscle', clwstrict=True)
        data = '\n'.join(['>%s\n' % str(sample_ids[i]) + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        stdout, stderr = muscle_cline(stdin=data)
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_vntrs = [str(aligned.seq) for aligned in alignment]

        return sample_ids, aligned_vntrs  # TODO sample_ids are not correctly sorted

    @staticmethod
    def _align_motifs_with_clustalo(sample_ids, labeled_vntrs, vid, score_matrix, output_dir):
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

    @staticmethod
    def _align_motifs_with_mafft(sample_ids, labeled_vntrs, vid, score_matrix, output_dir, preserve_order=False):
        # Check if MAFFT is installed
        if not shutil.which("mafft"):
            raise ValueError("MAFFT is not installed. Please install MAFFT and try again.")

        aln_input = f"{output_dir}/alignment_input.fa"
        aln_output = f"{output_dir}/alignment_output.fa"
        if vid is not None:
            aln_input = f"{output_dir}/{vid}_alignment_input.fa"
            aln_output = f"{output_dir}/{vid}_alignment_output.fa"

        # Writing alignment input
        data = '\n'.join(['>%s\n' % sample_ids[i] + labeled_vntrs[i] for i in range(len(labeled_vntrs))])
        with open(aln_input, "w") as f:
            f.write(data)

        def write_score_matrix_for_mafft(score_matrix, output_dir):
            score_matrix_file = "{}/matrixfile.txt".format(output_dir)
            of = open(score_matrix_file, "w")
            for symbol_1 in score_matrix:
                if symbol_1.startswith("gap"):
                    continue
                hex_symbol_1 = hex(ord(symbol_1))
                for symbol_2 in score_matrix[symbol_1]:
                    if symbol_2.startswith("gap"):
                        continue
                    hex_symbol_2 = hex(ord(symbol_2))
                    if symbol_2 == symbol_1:
                        # match
                        of.write("{} {} {}\n".format(hex_symbol_1, hex_symbol_2, score_matrix[symbol_1][symbol_2]))
                    else:
                        # mismatch
                        of.write("{} {} {}\n".format(hex_symbol_1, hex_symbol_2, score_matrix[symbol_1][symbol_2]))
            of.close()
            return score_matrix_file

        # TODO call mafft using pysam wrapper (
        # mafft_exe = "/usr/bin/mafft"
        # mafft_cline = MafftCommandline(mafft_exe, input=temp_input_name)

        if score_matrix is not None:
            score_matrix_file = write_score_matrix_for_mafft(score_matrix, output_dir)
            if preserve_order:
                os.system("mafft --quiet --auto "
                          "--textmatrix {} "
                          "--op {} --ep {} "
                          "{} > {}".format(score_matrix_file,
                                           score_matrix['gap_open'],
                                           score_matrix['gap_extension'],
                                           aln_input,
                                           aln_output))
            else:
                os.system("mafft --quiet --auto --reorder "
                # os.system("mafft --localpair --maxiterate 1000 --reorder "
                          "--textmatrix {} --op {} --ep {} "
                          "{} > {}".format(score_matrix_file,
                                           score_matrix['gap_open'],
                                           score_matrix['gap_extension'],
                                           aln_input,
                                           aln_output))
            os.remove(score_matrix_file)
        else:
            print("Score matrix file is not given. Default scoring parameters are used (not recommended).")
            if preserve_order:
                os.system("mafft --quiet --text --auto {} > {}".format(aln_input, aln_output))
            else:
                os.system("mafft --quiet --text --auto --reorder {} > {}".format(aln_input, aln_output))

        if not os.path.exists(aln_output):
            # os.system("mafft --text --auto {} > {}".format(aln_input, aln_output))  # to print out error
            raise FileNotFoundError("Error: motif alignment was not performed correctly.")

        alinged_sample_ids, aligned_vntrs = MotifAligner.load_aligned_trs(aln_output)

        return alinged_sample_ids, aligned_vntrs

    @staticmethod
    def load_aligned_trs(aln_output):
        """
        Loads aligned VNTRs from a file.

        Parameters
        ----------
        aln_output : str
            Path to the file containing aligned VNTRs.

        Returns
        -------
        alinged_sample_ids : list
            List of sample IDs.
        aligned_vntrs : list
            List of aligned VNTRs.
        """

        aligned_trs = []
        alinged_sample_ids = []
        tr_seq = ""
        with open(aln_output, "r") as f:
            for line in f:
                if line.startswith(">"):
                    alinged_sample_ids.append(line.strip()[1:])
                    if len(tr_seq) > 0:
                        aligned_trs.append(tr_seq)
                    tr_seq = ""
                else:
                    tr_seq += line.strip()
        if len(tr_seq) > 0:
            aligned_trs.append(tr_seq)
        if len(aligned_trs) == 0:
            raise ValueError(f"No aligned VNTRs in {aln_output}")
        return alinged_sample_ids, aligned_trs
