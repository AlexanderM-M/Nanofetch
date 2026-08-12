"""Dependency-free SVG coverage plots for NanoFetch extractions."""

import html
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

import pysam

from . import __version__
from .assemblies import resolve_contig
from .errors import NanoFetchError
from .models import GeneInterval, Region


@dataclass(frozen=True)
class CoverageTrack:
    region: Region
    primary: Tuple[float, ...]
    supplementary: Tuple[float, ...]
    secondary: Tuple[float, ...]

    @property
    def maximum(self) -> float:
        return max(
            (primary + supplementary + secondary
             for primary, supplementary, secondary in zip(
                 self.primary, self.supplementary, self.secondary
             )),
            default=0.0,
        )


@dataclass(frozen=True)
class CoverageSummary:
    symbol: str
    source: str
    assembly: str
    annotation: str
    gene_regions: Tuple[Region, ...]
    tracks: Tuple[CoverageTrack, ...]
    primary_alignments: int
    supplementary_alignments: int
    secondary_alignments: int
    aligned_bases: int

    @property
    def alignments(self) -> int:
        return (
            self.primary_alignments
            + self.supplementary_alignments
            + self.secondary_alignments
        )

    @property
    def span(self) -> int:
        return sum(track.region.end - track.region.start for track in self.tracks)

    @property
    def mean_depth(self) -> float:
        return self.aligned_bases / self.span if self.span else 0.0

    @property
    def maximum_depth(self) -> float:
        return max((track.maximum for track in self.tracks), default=0.0)


def _bin_bounds(region: Region, count: int, index: int) -> Tuple[int, int]:
    span = region.end - region.start
    return (
        region.start + index * span // count,
        region.start + (index + 1) * span // count,
    )


def _add_block(values: List[int], region: Region, start: int, end: int) -> int:
    """Add aligned bases from one block to bins and return the clipped total."""
    overlap_start = max(start, region.start)
    overlap_end = min(end, region.end)
    if overlap_start >= overlap_end:
        return 0
    count = len(values)
    span = region.end - region.start
    first = min(count - 1, (overlap_start - region.start) * count // span)
    last = min(count - 1, (overlap_end - 1 - region.start) * count // span)
    for index in range(first, last + 1):
        bin_start, bin_end = _bin_bounds(region, count, index)
        values[index] += max(0, min(overlap_end, bin_end) - max(overlap_start, bin_start))
    return overlap_end - overlap_start


def summarize_coverage(
    bam: Path,
    symbol: str,
    intervals: Sequence[GeneInterval],
    regions: Sequence[Region],
    assembly: str,
    annotation: str,
    source: str,
    bins: int = 240,
) -> CoverageSummary:
    """Calculate mean aligned-base depth in fixed-width bins from an output BAM."""
    if bins < 10:
        raise NanoFetchError("Coverage plots require at least 10 bins.")
    arrays = []
    for region in regions:
        count = min(bins, region.end - region.start)
        arrays.append({
            "primary": [0] * count,
            "supplementary": [0] * count,
            "secondary": [0] * count,
        })

    counts = {"primary": 0, "supplementary": 0, "secondary": 0}
    aligned_bases = 0
    with pysam.AlignmentFile(str(bam), "rb") as alignments:
        header_contigs = dict(zip(alignments.references, alignments.lengths))
        gene_regions = tuple(
            Region(
                resolve_contig(interval.contig, header_contigs),
                interval.start - 1,
                interval.end,
            )
            for interval in intervals
        )
        for alignment in alignments.fetch(until_eof=True):
            if alignment.is_unmapped:
                continue
            if alignment.is_supplementary:
                category = "supplementary"
            elif alignment.is_secondary:
                category = "secondary"
            else:
                category = "primary"
            counts[category] += 1
            contig = alignments.get_reference_name(alignment.reference_id)
            for block_start, block_end in alignment.get_blocks():
                for region, values in zip(regions, arrays):
                    if region.contig == contig:
                        aligned_bases += _add_block(
                            values[category], region, block_start, block_end
                        )

    tracks = []
    for region, values in zip(regions, arrays):
        count = len(values["primary"])
        widths = [
            _bin_bounds(region, count, index)[1]
            - _bin_bounds(region, count, index)[0]
            for index in range(count)
        ]
        tracks.append(CoverageTrack(
            region=region,
            primary=tuple(value / width for value, width in zip(values["primary"], widths)),
            supplementary=tuple(
                value / width for value, width in zip(values["supplementary"], widths)
            ),
            secondary=tuple(value / width for value, width in zip(values["secondary"], widths)),
        ))
    return CoverageSummary(
        symbol=symbol,
        source=source,
        assembly=assembly,
        annotation=annotation,
        gene_regions=gene_regions,
        tracks=tuple(tracks),
        primary_alignments=counts["primary"],
        supplementary_alignments=counts["supplementary"],
        secondary_alignments=counts["secondary"],
        aligned_bases=aligned_bases,
    )


def _format_number(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:.0f}"


def _format_coordinate(value: int) -> str:
    return f"{value / 1_000_000:.2f} Mb" if value >= 1_000_000 else f"{value:,} bp"


def _nice_ceiling(value: float) -> float:
    if value <= 1:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(value))
    fraction = value / magnitude
    step = 2 if fraction <= 2 else 5 if fraction <= 5 else 10
    return step * magnitude


def _area_path(lower: Sequence[float], upper: Sequence[float], x: float, y: float,
               width: float, height: float, scale: float) -> str:
    count = len(upper)
    points = []
    for index, value in enumerate(upper):
        left = x + width * index / count
        right = x + width * (index + 1) / count
        top = y + height - value / scale * height
        points.extend(((left, top), (right, top)))
    for index in range(count - 1, -1, -1):
        left = x + width * index / count
        right = x + width * (index + 1) / count
        bottom = y + height - lower[index] / scale * height
        points.extend(((right, bottom), (left, bottom)))
    if not points:
        return ""
    return "M " + " L ".join(f"{px:.2f} {py:.2f}" for px, py in points) + " Z"


def render_coverage_svg(summary: CoverageSummary) -> str:
    """Render a self-contained, accessible SVG document."""
    width = 1200
    left = 92
    plot_width = 1030
    plot_height = 190
    track_step = 270
    first_track_y = 220
    height = first_track_y + len(summary.tracks) * track_step + 62
    esc = lambda value: html.escape(str(value), quote=True)
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        'role="img" aria-labelledby="title description">',
        f'<title id="title">NanoFetch coverage for {esc(summary.symbol)}</title>',
        f'<desc id="description">Binned aligned-base depth across the extracted '
        f'region for {esc(summary.symbol)}.</desc>',
        '<style>',
        'text{font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif;fill:#0f172a}',
        '.muted{fill:#64748b}.small{font-size:13px}.label{font-size:14px;font-weight:600}',
        '.grid{stroke:#dbe3ee;stroke-width:1}.axis{stroke:#94a3b8;stroke-width:1.2}',
        '.primary{fill:#2563eb}.supplementary{fill:#f59e0b}.secondary{fill:#db2777}',
        '.gene{fill:#0f766e}.gene-zone{fill:#ccfbf1;opacity:.42}',
        '</style>',
        '<rect width="100%" height="100%" fill="#ffffff" rx="18"/>',
        f'<text x="{left}" y="52" font-size="26" font-weight="750">NanoFetch coverage</text>',
        f'<text x="{left}" y="82" font-size="19" font-weight="650">{esc(summary.symbol)}</text>',
        f'<text x="{left + 95}" y="82" class="muted" font-size="15">'
        f'{esc(summary.assembly)} · {esc(summary.annotation)} · {esc(summary.source)}</text>',
    ]

    cards = (
        ("ALIGNMENTS", f"{summary.alignments:,}"),
        ("MEAN DEPTH", f"{summary.mean_depth:.1f}×"),
        ("MAX BIN", f"{summary.maximum_depth:.1f}×"),
        ("EXTRACTED SPAN", _format_number(summary.span) + " bp"),
    )
    for index, (label, value) in enumerate(cards):
        x = left + index * 202
        lines.extend((
            f'<rect x="{x}" y="108" width="184" height="72" rx="10" fill="#f8fafc" '
            'stroke="#e2e8f0"/>',
            f'<text x="{x + 14}" y="132" class="muted small" font-weight="650">{label}</text>',
            f'<text x="{x + 14}" y="163" font-size="22" font-weight="720">{esc(value)}</text>',
        ))

    legend_x = 925
    for index, (label, color, count) in enumerate((
        ("Primary", "primary", summary.primary_alignments),
        ("Supplementary", "supplementary", summary.supplementary_alignments),
        ("Secondary", "secondary", summary.secondary_alignments),
    )):
        y = 115 + index * 24
        lines.append(f'<rect x="{legend_x}" y="{y}" width="12" height="12" rx="2" class="{color}"/>')
        lines.append(
            f'<text x="{legend_x + 19}" y="{y + 11}" class="small">'
            f'{label} {count:,}</text>'
        )

    for track_index, track in enumerate(summary.tracks):
        y = first_track_y + track_index * track_step
        region = track.region
        scale = _nice_ceiling(track.maximum)
        lines.append(
            f'<text x="{left}" y="{y - 16}" class="label">{esc(region.contig)}:'
            f'{region.start + 1:,}–{region.end:,}</text>'
        )
        for gene in summary.gene_regions:
            if gene.contig != region.contig:
                continue
            clipped_start = max(gene.start, region.start)
            clipped_end = min(gene.end, region.end)
            if clipped_start < clipped_end:
                gx = left + (clipped_start - region.start) / (region.end - region.start) * plot_width
                gw = (clipped_end - clipped_start) / (region.end - region.start) * plot_width
                lines.append(
                    f'<rect x="{gx:.2f}" y="{y}" width="{gw:.2f}" height="{plot_height}" '
                    'class="gene-zone"/>'
                )
        for fraction in (0, 0.5, 1):
            gy = y + plot_height - fraction * plot_height
            value = fraction * scale
            lines.append(
                f'<line x1="{left}" y1="{gy:.2f}" x2="{left + plot_width}" '
                f'y2="{gy:.2f}" class="grid"/>'
            )
            lines.append(
                f'<text x="{left - 12}" y="{gy + 5:.2f}" text-anchor="end" '
                f'class="muted small">{value:g}×</text>'
            )

        primary = track.primary
        primary_plus_supplementary = tuple(
            a + b for a, b in zip(track.primary, track.supplementary)
        )
        total = tuple(
            a + b + c for a, b, c in zip(
                track.primary, track.supplementary, track.secondary
            )
        )
        zero = (0.0,) * len(primary)
        for css_class, lower, upper in (
            ("primary", zero, primary),
            ("supplementary", primary, primary_plus_supplementary),
            ("secondary", primary_plus_supplementary, total),
        ):
            path = _area_path(lower, upper, left, y, plot_width, plot_height, scale)
            lines.append(f'<path d="{path}" class="{css_class}" opacity="0.82"/>')
        lines.append(
            f'<line x1="{left}" y1="{y + plot_height}" x2="{left + plot_width}" '
            f'y2="{y + plot_height}" class="axis"/>'
        )
        for fraction in (0, 0.5, 1):
            x = left + fraction * plot_width
            coordinate = round(region.start + fraction * (region.end - region.start))
            lines.append(
                f'<text x="{x:.2f}" y="{y + plot_height + 24}" text-anchor="middle" '
                f'class="muted small">{_format_coordinate(coordinate)}</text>'
            )

        lane_y = y + plot_height + 49
        lines.append(
            f'<line x1="{left}" y1="{lane_y}" x2="{left + plot_width}" '
            f'y2="{lane_y}" class="grid"/>'
        )
        for gene in summary.gene_regions:
            if gene.contig != region.contig:
                continue
            clipped_start = max(gene.start, region.start)
            clipped_end = min(gene.end, region.end)
            if clipped_start < clipped_end:
                gx = left + (clipped_start - region.start) / (region.end - region.start) * plot_width
                gw = max(2.0, (clipped_end - clipped_start) / (region.end - region.start) * plot_width)
                lines.append(
                    f'<rect x="{gx:.2f}" y="{lane_y - 7}" width="{gw:.2f}" height="14" '
                    'rx="4" class="gene"/>'
                )
                lines.append(
                    f'<text x="{gx + gw / 2:.2f}" y="{lane_y - 13}" text-anchor="middle" '
                    f'class="label">{esc(summary.symbol)}</text>'
                )

    lines.append(
        f'<text x="{left}" y="{height - 24}" class="muted small">'
        f'Binned aligned-base depth · generated by NanoFetch {esc(__version__)}</text>'
    )
    lines.append('</svg>')
    return "\n".join(lines) + "\n"


def write_coverage_svg(path: Path, summary: CoverageSummary, force: bool = False) -> None:
    """Atomically write a coverage SVG without replacing files by default."""
    if path.suffix.lower() != ".svg":
        raise NanoFetchError("--plot output must use the .svg extension.")
    if path.exists() and not force:
        raise NanoFetchError(f"Plot already exists: {path}. Use --force to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(render_coverage_svg(summary))
        os.replace(str(temporary), str(path))
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
