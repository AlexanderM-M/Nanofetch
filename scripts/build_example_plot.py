#!/usr/bin/env python3
"""Regenerate the SVG coverage example embedded in the README."""

from pathlib import Path

from nanofetch.models import Region
from nanofetch.plot import CoverageSummary, CoverageTrack, write_coverage_svg


def example_summary() -> CoverageSummary:
    bins = 120
    primary = tuple(
        max(0.0, 29.0 - abs(index - 60) * 0.42)
        + ((index * 7) % 6) * 0.35
        for index in range(bins)
    )
    supplementary = tuple(
        2.4 + (index % 4) * 0.25 if 30 <= index <= 93 and index % 9 < 4 else 0.0
        for index in range(bins)
    )
    secondary = tuple(
        1.8 if 43 <= index <= 84 and index % 13 < 3 else 0.0
        for index in range(bins)
    )
    region = Region("chr7", 54_018_819, 56_211_628)
    span = region.end - region.start
    return CoverageSummary(
        symbol="EGFR",
        source="example-tumor.bam",
        assembly="GRCh38",
        annotation="GENCODE 50 (GRCh38.p14)",
        gene_regions=(Region("chr7", 55_018_819, 55_211_628),),
        tracks=(CoverageTrack(region, primary, supplementary, secondary),),
        primary_alignments=18_432,
        supplementary_alignments=614,
        secondary_alignments=207,
        aligned_bases=round(14.7 * span),
    )


def main() -> None:
    destination = Path(__file__).resolve().parents[1] / "docs" / "example-coverage.svg"
    write_coverage_svg(destination, example_summary(), force=True)
    print(f"Wrote {destination}")


if __name__ == "__main__":
    main()
