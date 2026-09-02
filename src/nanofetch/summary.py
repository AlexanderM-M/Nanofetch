"""Compact per-gene coverage summaries."""

import csv
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence, Tuple

import pysam

from .errors import NanoFetchError
from .extract import iter_region_alignments, padded_regions
from .models import GeneInterval, Region


@dataclass(frozen=True)
class GeneSummary:
    symbol: str
    regions: Tuple[Region, ...]
    bases: int
    alignments: int
    aligned_bases: int
    covered_1x: int
    covered_10x: int
    covered_30x: int
    mapping_quality_sum: int
    mapping_quality_count: int

    @property
    def mean_depth(self) -> float:
        return self.aligned_bases / self.bases if self.bases else 0.0

    def covered_percent(self, count: int) -> float:
        return 100 * count / self.bases if self.bases else 0.0

    @property
    def mean_mapping_quality(self):
        if not self.mapping_quality_count:
            return None
        return self.mapping_quality_sum / self.mapping_quality_count


def summarize_gene(
    source: pysam.AlignmentFile,
    symbol: str,
    intervals: Sequence[GeneInterval],
    header_contigs: Mapping[str, int],
    include_supplementary: bool = False,
    include_secondary: bool = False,
) -> GeneSummary:
    """Calculate gene-body coverage using the extraction's record filters."""
    regions = padded_regions(intervals, header_contigs, 0)

    def selected(alignment: pysam.AlignedSegment) -> bool:
        return not (
            (alignment.is_secondary and not include_secondary)
            or (alignment.is_supplementary and not include_supplementary)
        )

    alignments = 0
    mapping_quality_sum = 0
    mapping_quality_count = 0
    for alignment in iter_region_alignments(source, regions):
        if not selected(alignment):
            continue
        alignments += 1
        if alignment.mapping_quality != 255:
            mapping_quality_sum += alignment.mapping_quality
            mapping_quality_count += 1

    bases = 0
    aligned_bases = 0
    covered_1x = 0
    covered_10x = 0
    covered_30x = 0
    for region in regions:
        bases += region.end - region.start
        coverage = source.count_coverage(
            region.contig,
            region.start,
            region.end,
            quality_threshold=0,
            read_callback=selected,
        )
        for depth in map(sum, zip(*coverage)):
            aligned_bases += depth
            covered_1x += depth >= 1
            covered_10x += depth >= 10
            covered_30x += depth >= 30

    return GeneSummary(
        symbol=symbol,
        regions=regions,
        bases=bases,
        alignments=alignments,
        aligned_bases=aligned_bases,
        covered_1x=covered_1x,
        covered_10x=covered_10x,
        covered_30x=covered_30x,
        mapping_quality_sum=mapping_quality_sum,
        mapping_quality_count=mapping_quality_count,
    )


def write_summary(
    path: Path, summaries: Sequence[GeneSummary], force: bool = False
) -> None:
    """Atomically write a compact TSV summary."""
    if path.exists() and not force:
        raise NanoFetchError(
            f"Summary already exists: {path}. Use --force to replace it."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "gene",
                    "regions",
                    "gene_bases",
                    "alignments",
                    "mean_depth",
                    "covered_1x_pct",
                    "covered_10x_pct",
                    "covered_30x_pct",
                    "mean_mapq",
                ]
            )
            for summary in summaries:
                regions = ",".join(
                    f"{region.contig}:{region.start + 1}-{region.end}"
                    for region in summary.regions
                )
                mean_mapq = summary.mean_mapping_quality
                writer.writerow(
                    [
                        summary.symbol,
                        regions,
                        summary.bases,
                        summary.alignments,
                        f"{summary.mean_depth:.3f}",
                        f"{summary.covered_percent(summary.covered_1x):.3f}",
                        f"{summary.covered_percent(summary.covered_10x):.3f}",
                        f"{summary.covered_percent(summary.covered_30x):.3f}",
                        "." if mean_mapq is None else f"{mean_mapq:.2f}",
                    ]
                )
        os.replace(str(temporary), str(path))
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
