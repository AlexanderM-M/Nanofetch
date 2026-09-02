import os
import tempfile
from itertools import chain
from pathlib import Path
from typing import Iterable, Iterator, List, Mapping, Sequence, Tuple

import pysam

from . import __version__
from .assemblies import resolve_contig
from .errors import NanoFetchError
from .models import ExtractionResult, GeneInterval, Region


def padded_regions(
    intervals: Sequence[GeneInterval],
    header_contigs: Mapping[str, int],
    padding: int,
) -> Tuple[Region, ...]:
    if padding < 0:
        raise NanoFetchError("--padding must be zero or greater.")
    raw: List[Region] = []
    for interval in intervals:
        contig = resolve_contig(interval.contig, header_contigs)
        length = header_contigs[contig]
        start = max(0, interval.start - 1 - padding)
        end = min(length, interval.end + padding)
        raw.append(Region(contig, start, end))

    return merge_regions(raw, header_contigs)


def merge_regions(
    regions: Iterable[Region], header_contigs: Mapping[str, int]
) -> Tuple[Region, ...]:
    """Return a coordinate-sorted union of regions."""
    order = {contig: index for index, contig in enumerate(header_contigs)}
    raw = sorted(
        regions, key=lambda region: (order[region.contig], region.start, region.end)
    )
    merged: List[Region] = []
    for region in raw:
        if (
            merged
            and merged[-1].contig == region.contig
            and region.start <= merged[-1].end
        ):
            previous = merged[-1]
            merged[-1] = Region(
                previous.contig, previous.start, max(previous.end, region.end)
            )
        else:
            merged.append(region)
    return tuple(merged)


def combined_regions(
    groups: Iterable[Sequence[Region]], header_contigs: Mapping[str, int]
) -> Tuple[Region, ...]:
    """Return the union of multiple genes' extraction regions."""
    return merge_regions(chain.from_iterable(groups), header_contigs)


def iter_region_alignments(
    source: pysam.AlignmentFile, regions: Sequence[Region]
) -> Iterator[pysam.AlignedSegment]:
    """Yield coordinate-sorted records once across a sorted region union."""
    previous = None
    for region in regions:
        for alignment in source.fetch(region.contig, region.start, region.end):
            end = alignment.reference_end
            if (
                previous is not None
                and previous.contig == region.contig
                and end is not None
                and alignment.reference_start < previous.end
                and end > previous.start
            ):
                continue
            yield alignment
        previous = region


def validate_input(source: pysam.AlignmentFile, path: Path) -> None:
    if not (source.is_bam or source.is_cram):
        raise NanoFetchError(
            f"Input must be BAM or CRAM; {path} is not a supported alignment file."
        )
    try:
        source.check_index()
    except (ValueError, OSError) as error:
        raise NanoFetchError(
            f"Input alignment is not indexed: {path}. Create an index with "
            f"'samtools index {path}' or 'pysam.index(\"{path}\")'."
        ) from error
    sort_order = source.header.to_dict().get("HD", {}).get("SO")
    if sort_order and sort_order != "coordinate":
        raise NanoFetchError(
            f"Input alignment declares sort order {sort_order!r}; "
            "coordinate sorting is required."
        )


def header_with_program(source: pysam.AlignmentFile) -> dict:
    header = source.header.to_dict()
    programs = header.setdefault("PG", [])
    used = {record.get("ID") for record in programs}
    program_id = "nanofetch"
    suffix = 1
    while program_id in used:
        suffix += 1
        program_id = f"nanofetch.{suffix}"
    programs.append({"ID": program_id, "PN": "nanofetch", "VN": __version__})
    return header


def output_name(symbol: str) -> str:
    safe = "".join(char if char.isalnum() or char in "._-" else "_" for char in symbol)
    return f"{safe}.bam"


def extract_gene(
    source: pysam.AlignmentFile,
    symbol: str,
    intervals: Sequence[GeneInterval],
    output_dir: Path,
    padding: int = 1_000_000,
    include_supplementary: bool = False,
    include_secondary: bool = False,
    threads: int = 1,
    index_output: bool = False,
    force: bool = False,
) -> ExtractionResult:
    header_contigs = dict(zip(source.references, source.lengths))
    regions = padded_regions(intervals, header_contigs, padding)
    output = output_dir / output_name(symbol)
    return extract_regions(
        source,
        symbol,
        regions,
        output,
        include_supplementary=include_supplementary,
        include_secondary=include_secondary,
        threads=threads,
        index_output=index_output,
        force=force,
    )


def extract_regions(
    source: pysam.AlignmentFile,
    symbol: str,
    regions: Sequence[Region],
    output: Path,
    include_supplementary: bool = False,
    include_secondary: bool = False,
    threads: int = 1,
    index_output: bool = False,
    force: bool = False,
) -> ExtractionResult:
    """Extract a sorted region union to one BAM."""
    if threads < 1:
        raise NanoFetchError("--threads must be at least 1.")
    if output.exists() and not force:
        raise NanoFetchError(
            f"Output already exists: {output}. Use --force to replace it."
        )
    output.parent.mkdir(parents=True, exist_ok=True)

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=str(output.parent)
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    count = 0
    try:
        with pysam.AlignmentFile(
            str(temporary), "wb", header=header_with_program(source), threads=threads
        ) as destination:
            for alignment in iter_region_alignments(source, regions):
                if alignment.is_secondary and not include_secondary:
                    continue
                if alignment.is_supplementary and not include_supplementary:
                    continue
                destination.write(alignment)
                count += 1
        os.replace(str(temporary), str(output))
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    # A replaced BAM must never retain an index for its previous contents.
    for suffix in (".bai", ".csi"):
        Path(str(output) + suffix).unlink(missing_ok=True)
    index_path = None
    if index_output:
        try:
            args = ["-@", str(threads), str(output)] if threads > 1 else [str(output)]
            pysam.index(*args)
            index_path = Path(str(output) + ".bai")
        except Exception as error:
            raise NanoFetchError(
                f"Created {output}, but indexing failed: {error}"
            ) from error
    return ExtractionResult(symbol, output, regions, count, index_path)
