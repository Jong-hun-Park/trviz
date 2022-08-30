import pytest

from trviz.decomposer import TandemRepeatDecomposer


@pytest.fixture(scope="session")
def tr_decomposer_dp():
    return TandemRepeatDecomposer(mode="DP")


@pytest.mark.parametrize(
    "sequence, motifs, expected",
    [
        (
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "AAAAAA",
                ["AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA"]
        ),
        (
                "ACTGACTGACTG",
                "ACTG",
                ["ACTG", "ACTG", "ACTG"]
        ),
        (
                "AACCTTTTCTAACCTTTTCT",
                "AACCTTTTCT",
                ["AACCTTTTCT", "AACCTTTTCT"]
        ),
        (
                "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG",
                "CGG",
                ["CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                 "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG"]
        ),
    ]
)
def test_decompose_dp_perfect_repeats_single_motif(tr_decomposer_dp, sequence, motifs, expected):
    decomposed_tr = tr_decomposer_dp.decompose(sequence, motifs)
    assert decomposed_tr == expected


@pytest.mark.parametrize(
    "sequence, motifs, expected",
    [
        (
                "AAAAAC"
                "AAAAAA"
                "AAAAAT"
                "AAAAAA"
                "TTAAAA",
                ["AAAAAA"],
                ["AAAAAC", "AAAAAA", "AAAAAT", "AAAAAA", "TTAAAA"]
        ),
        (
                "ACTG"
                "ACTT"
                "ACTG",
                ["ACTG"],
                ["ACTG", "ACTT", "ACTG"]
        ),
        (
                "AACCTTTTCT"
                "AACCTTGTCT",
                ["AACCTTTTCT"],
                ["AACCTTTTCT", "AACCTTGTCT"]
        ),
        (
                "CGCCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGT",
                ["CGG"],
                # ["CGC", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                #  "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGT"]
                ["CG", "CCGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                 "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGT"]
        ),
        (
                "AAATA"
                "AAATT"
                "AAATA"
                "AAAATA",
                ["AAATA"],
                ['AAATA', 'AAATT', 'AAATA', 'AAAATA']
        ),
    ]
)
def test_decompose_dp_imperfect_repeats_single_motif(tr_decomposer_dp, sequence, motifs, expected):
    decomposed_tr = tr_decomposer_dp.decompose(sequence, motifs)
    assert decomposed_tr == expected


@pytest.mark.parametrize(
    "sequence, motifs, expected",
    [
        (
                "AAAAAC"
                "AAAAAA"
                "AAAAAT"
                "AAAAAA"
                "TTAAAA",
                ["AAAAAA", "TTAAAA"],
                ["AAAAAC", "AAAAAA", "AAAAAT", "AAAAAA", "TTAAAA"]
        ),
        (
                "ACTG"
                "ACTT"
                "ACTG",
                ["ACTG", "ACTT"],
                ["ACTG", "ACTT", "ACTG"]
        ),
        (
                "AACCTTTTCT"
                "AACCTTGTCT"
                "AACCTTGTCT",
                ["AACCTTTTCT", "AACCTTGTCT"],
                ["AACCTTTTCT", "AACCTTGTCT", "AACCTTGTCT"]
        ),
        (
                "CGCCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGT",
                ["CGG", "CGC", "CGT"],
                ["CGC", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                 "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGT"]
                # Compare with the prev test
                # ["CG", "CCGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                #  "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGT"]
        ),
        (
                "AAATA"
                "AAATT"
                "AAATA"
                "AAAATA",
                ["AAATA", "AAATAA"],
                # ['AAATA', 'AAATT', 'AAATA', 'AAAATA']  # Compare with the prev test
                ['AAATA', 'AAATT', 'AAATAA', 'AAATA']
        ),
    ]
)
def test_decompose_dp_imperfect_repeats_multiple_motif(tr_decomposer_dp, sequence, motifs, expected):
    decomposed_tr = tr_decomposer_dp.decompose(sequence, motifs)
    assert decomposed_tr == expected


def test_decompose_dp_invalid_sequence(tr_decomposer_dp):
    with pytest.raises(ValueError):
        tr_decomposer_dp.decompose("NNNNNN", "ACTG")
