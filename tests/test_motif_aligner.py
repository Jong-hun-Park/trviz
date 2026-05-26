import pytest

from trviz.motif_aligner import MotifAligner


def test_get_motif_aligner_rejects_unknown_tool():
    aligner = MotifAligner()
    with pytest.raises(ValueError, match="Unknown alignment tool"):
        aligner._get_motif_aligner("not-a-real-tool")
