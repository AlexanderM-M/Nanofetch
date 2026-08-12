import os
import tempfile
from pathlib import Path
from typing import List, Mapping, Sequence, Tuple

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

    order = {contig: index for index, contig in enumerate(header_contigs)}
    raw.sort(key=lambda region: (order[region.contig], region.start, region.end))
    merged: List[Region] = []
    for region in raw:
        if merged and merged[-1].contig == region.contig and region.start <= merged[-1].end:
            previous = merged[-1]
            merged[-1] = Region(previous.contig, previous.start, max(previous.end, region.end))
        else:
            merged.append(region)
    return tuple(merged)


def validate_input(source: pysam.AlignmentFile, path: Path) -> None:
    if not source.is_bam:
        raise NanoFetchError(f"Input must be BAM; {path} is not a BAM file.")
    try:
        source.check_index()
    except (ValueError, OSError) as error:
        raise NanoFetchError(
            f"Input BAM is not indexed: {path}. Create an index with "
            f"'samtools index {path}' or 'pysam.index(\"{path}\")'."
        ) from error
    sort_order = source.header.to_dict().get("HD", {}).get("SO")
    if sort_order and sort_order != "coordinate":
        raise NanoFetchError(
            f"Input BAM declares sort order {sort_order!r}; coordinate sorting is required."
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
    if threads < 1:
        raise NanoFetchError("--threads must be at least 1.")
    header_contigs = dict(zip(source.references, source.lengths))
    regions = padded_regions(intervals, header_contigs, padding)
    output = output_dir / output_name(symbol)
    if output.exists() and not force:
        raise NanoFetchError(f"Output already exists: {output}. Use --force to replace it.")
    output_dir.mkdir(parents=True, exist_ok=True)

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=str(output_dir)
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    count = 0
    try:
        with pysam.AlignmentFile(
            str(temporary), "wb", header=header_with_program(source), threads=threads
        ) as destination:
            for region in regions:
                for alignment in source.fetch(region.contig, region.start, region.end):
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
            raise NanoFetchError(f"Created {output}, but indexing failed: {error}") from error
    return ExtractionResult(symbol, output, regions, count, index_path)
