from dataclasses import replace
from pathlib import Path
from runpy import run_path
from xml.etree import ElementTree

import pysam
import pytest

from nanofetch.annotations import resolve_gene
from nanofetch.cli import run
from nanofetch.errors import NanoFetchError
from nanofetch.extract import extract_gene
from nanofetch.plot import render_coverage_svg, summarize_coverage


EXAMPLE_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_example_plot.py"
example_summary = run_path(str(EXAMPLE_SCRIPT))["example_summary"]


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


def test_plot_explains_depth_and_bin_resolution():
    svg = render_coverage_svg(example_summary())
    assert "MEAN REGION DEPTH" in svg
    assert "Mean depth includes padding" in svg
    assert "maximum is a bin average" in svg
    assert "120 bins · ~18.3k bp/bin" in svg


def test_plot_wraps_long_metadata_above_summary_metrics():
    source = "sample_" + "a" * 240 + "<&>.bam"
    summary = replace(example_summary(), symbol="A_LONG_GENE_SYMBOL", source=source)
    root = ElementTree.fromstring(render_coverage_svg(summary))
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    headings = [
        node for node in root.findall("svg:text", namespace)
        if node.get("font-size") == "15"
    ]
    assert len(headings) > 1
    assert source in "".join(node.text for node in headings)
    metrics = root.find("svg:g[@id='summary-metrics']", namespace)
    divider = metrics.find("svg:line", namespace)
    assert max(float(node.get("y")) for node in headings) < float(divider.get("y1"))
