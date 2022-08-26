from collections import Counter
from typing import Dict

from trviz.utils import INDEX_TO_CHR, PRIVATE_MOTIF_LABEL
from trviz.utils import get_motif_counter


class MotifEncoder:

    def __init__(self, private_motif_threshold=0):
        self.private_motif_threshold = private_motif_threshold

    @staticmethod
    def _divide_motifs_into_normal_and_private(motif_counter, private_motif_threshold):
        """
        Givne a list of decomposed VNTRs, divide motifs into two groups: normal and private.
        If a motif occurred less than the private_motif_threshold, it is regarded as private motif.
        Otherwise, normal motifs.

        :param decomposed_vntrs
        :param private_motif_threshold
        :return: normal motifs, private motifs
        """

        normal_motifs = []
        private_motifs = []
        for motif, count in motif_counter.most_common():
            if count > private_motif_threshold:
                normal_motifs.append(motif)
            else:
                private_motifs.append(motif)

        return normal_motifs, private_motifs

    @staticmethod
    def find_minimum_private_motif_threshold(decomposed_vntrs):
        motif_counter = get_motif_counter(decomposed_vntrs)

        min_private_motif_threshold = 0
        for index, (motif, count) in enumerate(motif_counter.most_common()):
            print(index, motif, count)
            if index + 1 > len(INDEX_TO_CHR) - 1:  # starting with 33, last 126, skipping 4 symbols, 1 for private
                min_private_motif_threshold = count
                break

        return min_private_motif_threshold

    @staticmethod
    def write_motif_map(output_file, motif_to_alphabet, motif_counter):
        with open(output_file, "w") as f:
            for (motif, _) in motif_counter.most_common():
                f.write(f"{motif}\t{motif_to_alphabet[motif]}\t{motif_counter[motif]}\n")

    @staticmethod
    def _encode_decomposed_vntr(decomposed_vntrs, motif_to_symbol):
        labeled_trs = []
        for vntr in decomposed_vntrs:
            labeled_vntr = ""
            for motif in vntr:
                labeled_vntr += str(motif_to_symbol[motif])
            labeled_trs.append(labeled_vntr)
        return labeled_trs

    def encode(self, decomposed_vntrs, motif_map_file, auto=True):
        """

        :param decomposed_vntrs:
        :param private_motif_threshold:
        :param auto: if True, find the minimum threshold to encode everything using 126-33+1 ASCII characters.
        :return: labeled_vntrs, motif to alphabet (dictionary of the mapping)
        """

        def _index_to_char(index):
            """
            --anysymbol
            To use unusual characters (e.g., U as selenocysteine in protein sequence; i as inosine in nucleotide sequence),
            use the --anysymbol option:
            % mafft --anysymbol input > output
            It accepts any printable characters (U, O, #, $, %, etc.; 0x21-0x7e in the ASCII code),
            execpt for > (0x3e) and ( (0x28).
            # '= (60) < (61) > (62)' can not be used
            Unusual characters are scored as unknown (not considered in the calculation), unlike in the --text mode.
            """
            if index < 0 or index > len(INDEX_TO_CHR) - 1:
                raise ValueError(f"Index should range between 0 to {len(INDEX_TO_CHR) - 1}. Given : {index}")

            return INDEX_TO_CHR[index]

        if auto:
            self.private_motif_threshold = self.find_minimum_private_motif_threshold(decomposed_vntrs)
            print("private motif threshold: ", self.private_motif_threshold)

        motif_to_symbol: Dict = {}
        symbol_to_motif: Dict = {}
        motif_counter: Counter = get_motif_counter(decomposed_vntrs)

        # For private motifs, we use single letter to encode them.
        if self.private_motif_threshold > 0:
            normal_motifs, private_motifs = self._divide_motifs_into_normal_and_private(motif_counter,
                                                                                        self.private_motif_threshold)
            if len(normal_motifs) + 1 > len(INDEX_TO_CHR):
                print("Motif counter:", motif_counter)
                raise ValueError("Too many unique motifs. Can not encode properly: {} unique motifs".format(
                    len(normal_motifs) + len(private_motifs)))

            # Assign a code to all private motifs
            motif_to_symbol.update({motif: PRIVATE_MOTIF_LABEL for motif in sorted(private_motifs)})
            symbol_to_motif.update({PRIVATE_MOTIF_LABEL: motif for motif in sorted(private_motifs)})

            # For normal motifs
            motif_to_symbol.update(
                {motif: _index_to_char(index) for index, motif in enumerate(sorted(normal_motifs))})
            symbol_to_motif.update(
                {_index_to_char(index): motif for index, motif in enumerate(sorted(normal_motifs))})
        else:  # Use all distinct motif
            if (unique_motif_count := len(motif_counter)) > len(INDEX_TO_CHR):  # single symbol ascii
                raise ValueError(
                    "Too many unique motifs. Can not encode properly: {} unique motifs".format(unique_motif_count))

            motif_to_symbol.update(
                {motif: _index_to_char(index) for index, (motif, _) in enumerate(motif_counter.most_common())})
            symbol_to_motif.update(
                {_index_to_char(index): motif for index, (motif, _) in enumerate(motif_counter.most_common())})

        # Write motif encoding
        self.write_motif_map(motif_map_file, motif_to_symbol, motif_counter)

        # Encode TRs
        encoded_vntrs = self._encode_decomposed_vntr(decomposed_vntrs, motif_to_symbol)

        return encoded_vntrs

