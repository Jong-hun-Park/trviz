"""Tests for the matplotlib-idiomatic fig/ax return from the visualization API."""

from io import BytesIO

import matplotlib

matplotlib.use("Agg")  # headless backend for CI

import matplotlib.pyplot as plt
import pytest

from trviz.visualizer import TandemRepeatVisualizer


@pytest.fixture
def visualizer():
    return TandemRepeatVisualizer()


@pytest.fixture
def simple_inputs():
    """Two aligned alleles encoded as motif symbols."""
    aligned = ["AAAB", "AABB"]  # encoded post-alignment
    sample_ids = ["s1", "s2"]
    symbol_to_motif = {"A": "ACT", "B": "ACG"}
    return aligned, sample_ids, symbol_to_motif


def test_trplot_returns_fig_and_ax(visualizer, simple_inputs):
    aligned, sample_ids, symbol_to_motif = simple_inputs
    fig, ax = visualizer.trplot(
        aligned_labeled_repeats=aligned,
        sample_ids=sample_ids,
        symbol_to_motif=symbol_to_motif,
        sort_by_clustering=False,
        output_name=None,  # skip savefig
    )
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    plt.close(fig)


def test_trplot_post_hoc_styling_works(visualizer, simple_inputs):
    """The whole point: users can modify the returned axes/figure after the call."""
    aligned, sample_ids, symbol_to_motif = simple_inputs
    fig, ax = visualizer.trplot(
        aligned_labeled_repeats=aligned,
        sample_ids=sample_ids,
        symbol_to_motif=symbol_to_motif,
        sort_by_clustering=False,
        output_name=None,
    )
    # post-hoc adjustments — the user's main use case
    ax.set_title("Adjusted title")
    ax.tick_params(labelsize=14)
    assert ax.get_title() == "Adjusted title"

    # render to a buffer to confirm the figure is still alive
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=72)
    assert buf.tell() > 0  # bytes were written
    plt.close(fig)


def test_trplot_does_not_close_by_default(visualizer, simple_inputs):
    """Backward-incompatible behavior change: plt.close() is no longer automatic."""
    aligned, sample_ids, symbol_to_motif = simple_inputs
    fig, _ = visualizer.trplot(
        aligned_labeled_repeats=aligned,
        sample_ids=sample_ids,
        symbol_to_motif=symbol_to_motif,
        sort_by_clustering=False,
        output_name=None,
    )
    # If close were happening, fig.number would be missing from plt.get_fignums()
    assert fig.number in plt.get_fignums()
    plt.close(fig)


def test_trplot_closes_when_close_true(visualizer, simple_inputs):
    aligned, sample_ids, symbol_to_motif = simple_inputs
    fig, _ = visualizer.trplot(
        aligned_labeled_repeats=aligned,
        sample_ids=sample_ids,
        symbol_to_motif=symbol_to_motif,
        sort_by_clustering=False,
        output_name=None,
        close=True,
    )
    assert fig.number not in plt.get_fignums()


def test_trplot_with_user_provided_ax(visualizer, simple_inputs):
    """Users can inject their own ax for custom subplot layouts."""
    aligned, sample_ids, symbol_to_motif = simple_inputs
    user_fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    returned_fig, returned_ax = visualizer.trplot(
        aligned_labeled_repeats=aligned,
        sample_ids=sample_ids,
        symbol_to_motif=symbol_to_motif,
        sort_by_clustering=False,
        output_name=None,
        ax=axes[0],
    )
    assert returned_fig is user_fig
    assert returned_ax is axes[0]
    plt.close(user_fig)


def test_trplot_with_single_sequence(visualizer):
    """Regression: clustering is a no-op for N=1, must not crash on empty distance matrix."""
    fig, ax = visualizer.trplot(
        aligned_labeled_repeats=["AAAB"],
        sample_ids=["s1"],
        symbol_to_motif={"A": "ACT", "B": "ACG"},
        sort_by_clustering=True,  # the buggy path — was raising ValueError
        output_name=None,
    )
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    plt.close(fig)


def test_plot_motif_color_map_returns_fig_and_ax(visualizer):
    symbol_to_motif = {"A": "ACT", "B": "ACG"}
    motif_counter = {"ACT": 3, "ACG": 1}
    symbol_to_color = {"A": (0.9, 0.6, 0.0), "B": (0.3, 0.7, 0.9)}
    fig, ax = visualizer.plot_motif_color_map(
        symbol_to_motif,
        motif_counter,
        symbol_to_color,
        file_name=None,  # skip savefig
    )
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    # confirm post-hoc styling works
    ax.set_title("My motif map")
    assert ax.get_title() == "My motif map"
    plt.close(fig)
