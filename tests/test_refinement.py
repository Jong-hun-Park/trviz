import pytest
from trviz.decomposer import Decomposer


@pytest.mark.parametrize(
    "decomposed_trs, expected",
    [
        (
            [['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
             ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
             ['AACA', 'TAACA', 'AACA', 'AACA', 'TAACA']],
            [['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
             ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
             ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA']],
        ),
    ]
)
def test_refine(decomposed_trs, expected):
    refined_trs = Decomposer.refine(decomposed_trs)
    assert refined_trs == expected