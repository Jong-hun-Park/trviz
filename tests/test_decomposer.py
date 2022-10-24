import pytest

from trviz.decomposer import Decomposer


@pytest.fixture(scope="session")
def tr_decomposer_dp():
    return Decomposer(mode="DP")


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
                ["CGC", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                 "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGT"]
                # ["CG", "CCGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                #  "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGT"]
        ),
        (
                "AAATA"
                "AAATT"
                "AAATAA"
                "AAATA",
                ["AAATA"],
                ['AAATA', 'AAATT', 'AAATAA', 'AAATA']
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
        (
                "ACCCA"
                "ACCC"
                "ACCCA"
                "ACCCA",
                ["ACCCA"],
                ["ACCCA", "ACCC", "ACCCA", "ACCCA"],
        ),
        (
                "ACT"
                "ACT"
                "ACC"
                "ACT",
                ["ACT"],
                ["ACT", "ACT", "ACC", "ACT"],
        ),
        (
                "ACT"
                "ACT"
                "ACCG"
                "ACT",
                ["ACT"],
                ["ACT", "ACT", "ACCG", "ACT"],
        ),
        (
                "ACT"
                "ACT"
                "AC"
                "CG"
                "ACT",
                ["ACT", "AC", "CG"],
                ["ACT", "ACT", "AC", "CG", "ACT"],
        ),
        (
                "AATAA"
                "AATAAA"
                "AATAA",
                ["AATAA"],
                ["AATAA", "AATAAA", "AATAA"],
        ),
    ]
)
def test_decompose_dp_imperfect_repeats_multiple_motif(tr_decomposer_dp, sequence, motifs, expected):
    decomposed_tr = tr_decomposer_dp.decompose(sequence, motifs)
    assert decomposed_tr == expected


@pytest.mark.parametrize(
    "sequence, motifs, kwargs, expected",
    [
        (
                "ACGTTTACGTTTACGTTTACGTTT",
                ["ACGTTT"],
                {'match_score': 5, 'mismatch_score': -2, 'insertion_score': -3, 'deletion_score': -3},
                ["ACGTTT", "ACGTTT", "ACGTTT", "ACGTTT"],
        ),
        (
                "ACGTTTACGTTTACGTTTACGTTT",
                "ACGTTT",
                {'match_score': '5', 'mismatch_score': '3', 'insertion_score': '1', 'deletion_score': '1'},
                ValueError,
        ),
        (
                "ACGTTTACGTTTACGTTTACGTTT",
                "ACGTTT",
                {'mismatch': 3},  # Invalid. mismatch_score is the correct key
                KeyError,
        ),
        (
                "ACGTTTACGTTTACGTTTACGTTT",
                "ACGTTT",
                {'verbose': True},
                ["ACGTTT", "ACGTTT", "ACGTTT", "ACGTTT"],
        ),
        (
                "ACGTTTACGTTTACGTTTACGTTT",
                "ACGTTT",
                {'verbose': 'T'},
                ValueError,
        ),
    ]
)
def test_decompose_dp_arguments(tr_decomposer_dp, sequence, motifs, kwargs, expected):
    if type(expected) == type and issubclass(expected, Exception):
        with pytest.raises(expected):
            tr_decomposer_dp.decompose(sequence, motifs, **kwargs)
    else:
        decomposed_tr = tr_decomposer_dp.decompose(sequence, motifs, **kwargs)
        assert decomposed_tr == expected


def test_decompose_dp_invalid_sequence(tr_decomposer_dp):
    with pytest.raises(ValueError):
        tr_decomposer_dp.decompose("NNNNNN", "ACTG")


@pytest.mark.parametrize(
    "sequence, motifs, kwargs, expected",
    [
        (  # insertion tie-break
                "ACGTACGTTACGTAACGT",
                ["ACGT"],
                {'match_score': 2, 'mismatch_score': -1, 'insertion_score': -1, 'deletion_score': -1, 'verbose': True},
                ["ACGT", "ACGTT", "ACGTA", "ACGT"],
        ),
        (  # insertion tie-break
                "ACGTACGTTACGTAACGT",
                ["ACGT"],
                {'match_score': 2, 'mismatch_score': -1, 'insertion_score': -2, 'deletion_score': -2, 'verbose': True},
                ["ACGT", "ACGTT", "ACGTA", "ACGT"],
        ),
        (
                "ACCCAACCCACCCAACCCA",
                ["ACCCA"],
                {'match_score': 2, 'mismatch_score': -1, 'insertion_score': -1, 'deletion_score': -1, 'verbose': True},
                ["ACCCA", "ACCC", "ACCCA", "ACCCA"],
        ),
    ]
)
def test_decompose_dp_tie_break(tr_decomposer_dp, sequence, motifs, kwargs, expected):
    if type(expected) == type and issubclass(expected, Exception):
        with pytest.raises(expected):
            tr_decomposer_dp.decompose(sequence, motifs, **kwargs)
    else:
        decomposed_tr = tr_decomposer_dp.decompose(sequence, motifs, **kwargs)
        assert decomposed_tr == expected


