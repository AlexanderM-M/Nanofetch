from pathlib import Path

import pysam
import pytest

from nanofetch.annotations import resolve_gene
from nanofetch.cli import run
from nanofetch.errors import NanoFetchError
from nanofetch.extract import extract_gene
from nanofetch.plot import render_coverage_svg, summarize_coverage
from scripts.build_example_plot import example_summary


def test_coverage_summary_counts_alignment_classes(grch38_bam, tmp_path):
    with pysam.AlignmentFile(grch38_bam, "rb") as source:
        result = extract_gene(
            source,
            "EGFR",
            resolve_gene("EGFR", "grch38"),
            tmp_path / "output",
            padding=0,
            include_supplementary=True,
            include_secondary=True,
        )
    summary = summarize_coverage(
        result.output,
        "EGFR",
        resolve_gene("EGFR", "grch38"),
        result.regions,
        "GRCh38",
        "GENCODE 50 (GRCh38.p14)",
        "input.bam",
        bins=20,
    )
    assert summary.primary_alignments == 1
    assert summary.supplementary_alignments == 1
    assert summary.secondary_alignments == 1
    assert summary.aligned_bases == 300
    assert summary.maximum_depth > 0


def test_cli_writes_self_contained_svg(grch38_bam, tmp_path):
    plot = tmp_path / "plots" / "EGFR.coverage.svg"
    assert run([
        str(grch38_bam),
        "EGFR",
        "--padding", "0",
        "--output-dir", str(tmp_path / "output"),
        "--include-supplementary",
        "--include-secondary",
        "--plot", str(plot),
    ]) == 0
    svg = plot.read_text(encoding="utf-8")
    assert svg.startswith('<?xml version="1.0" encoding="UTF-8"?>')
    assert '<svg xmlns="http://www.w3.org/2000/svg"' in svg
    assert "NanoFetch coverage for EGFR" in svg
    assert "Primary 1" in svg
    assert "Supplementary 1" in svg
    assert "Secondary 1" in svg
    assert "GENCODE 50 (GRCh38.p14)" in svg


def test_plot_requires_one_gene_and_svg_extension(grch38_bam, tmp_path):
    with pytest.raises(NanoFetchError, match="exactly one"):
        run([
            str(grch38_bam), "EGFR", "MET", "--padding", "0",
            "--output-dir", str(tmp_path / "output"),
            "--plot", str(tmp_path / "plot.svg"),
        ])
    with pytest.raises(NanoFetchError, match=".svg extension"):
        run([
            str(grch38_bam), "EGFR", "--padding", "0",
            "--output-dir", str(tmp_path / "output"),
            "--plot", str(tmp_path / "plot.png"),
        ])


def test_documented_example_matches_renderer():
    example = Path(__file__).resolve().parents[1] / "docs" / "example-coverage.svg"
    assert example.read_text(encoding="utf-8") == render_coverage_svg(example_summary())
