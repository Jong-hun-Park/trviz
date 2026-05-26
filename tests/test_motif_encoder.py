"""Tests for trviz.motif_encoder.MotifEncoder."""

from trviz.motif_encoder import MotifEncoder
from trviz.utils import PRIVATE_MOTIF_LABEL


def test_private_motifs_attribute_lists_all_collapsed_motifs(tmp_path):
    """When multiple distinct motifs share PRIVATE_MOTIF_LABEL, MotifEncoder
    must expose all of them via .private_motifs so the legend can label the
    private color correctly. symbol_to_motif itself only holds one arbitrary
    private motif (intentional: many → one in the encoded space)."""
    decomposed_trs = [
        ["ACT"] * 4 + ["AAA"],
        ["ACT"] * 4 + ["CCC"],
        ["ACT"] * 4 + ["GGG"],
        ["ACT"] * 5,
    ]
    encoder = MotifEncoder(private_motif_threshold=2)
    encoder.encode(decomposed_trs, motif_map_file=str(tmp_path / "map.txt"), auto=False)

    # All 3 single-occurrence motifs are private (count < threshold of 2)
    assert set(encoder.private_motifs) == {"AAA", "CCC", "GGG"}

    # symbol_to_motif still uses the lossy many-to-one shape
    assert encoder.symbol_to_motif[PRIVATE_MOTIF_LABEL] in {"AAA", "CCC", "GGG"}

    # All 3 private motifs map to the same symbol (the encoding's intent)
    private_symbols = {encoder.motif_to_symbol[m] for m in ("AAA", "CCC", "GGG")}
    assert private_symbols == {PRIVATE_MOTIF_LABEL}


def test_private_motifs_empty_when_no_private_threshold(tmp_path):
    """With private_motif_threshold=0, no motifs are private and the list is empty."""
    decomposed_trs = [["ACT", "ACT"], ["ACT", "ACG"]]
    encoder = MotifEncoder(private_motif_threshold=0)
    encoder.encode(decomposed_trs, motif_map_file=str(tmp_path / "map.txt"), auto=False)
    assert encoder.private_motifs == []
