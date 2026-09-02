import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path
from typing import List, Sequence, Tuple

import pysam

from . import __version__
from .annotations import annotation_metadata, resolve_gene, unique_symbols
from .assemblies import ASSEMBLY_LABELS, choose_assembly
from .bed import write_bed
from .errors import NanoFetchError
from .extract import (
    combined_regions,
    extract_gene,
    extract_regions,
    output_name,
    padded_regions,
    validate_input,
)
from .panels import available_panels, load_panel, panel_descriptions
from .plot import summarize_coverage, write_coverage_svg
from .summary import summarize_gene, write_summary


def nonnegative_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if number < 0:
        raise argparse.ArgumentTypeError("must be zero or greater")
    return number


def positive_integer(value: str) -> int:
    number = nonnegative_integer(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="nanofetch",
        description="Ridiculously easy gene-region alignment extraction.",
    )
    result.add_argument(
        "input", nargs="?", type=Path, help="coordinate-sorted, indexed BAM or CRAM"
    )
    result.add_argument("gene", nargs="*", help="one or more gene symbols")
    result.add_argument(
        "--genes", nargs="+", default=[], metavar="GENE", help="additional gene symbols"
    )
    result.add_argument(
        "--gene-file",
        action="append",
        default=[],
        type=Path,
        metavar="FILE",
        help="gene symbols or compact families, with # comments",
    )
    result.add_argument(
        "--panel", action="append", default=[], metavar="NAME", help="built-in panel"
    )
    result.add_argument(
        "--list-panels", action="store_true", help="list built-in panels and exit"
    )
    result.add_argument(
        "--genome",
        default="auto",
        metavar="BUILD",
        help="auto, grch37/hg19, grch38/hg38, or t2t/hs1 (default: auto)",
    )
    result.add_argument(
        "--reference", type=Path, metavar="FASTA", help="reference FASTA for CRAM input"
    )
    result.add_argument(
        "--padding", type=nonnegative_integer, default=1_000_000, metavar="BP"
    )
    result.add_argument("--include-supplementary", action="store_true")
    result.add_argument("--include-secondary", action="store_true")
    result.add_argument("--threads", type=positive_integer, default=1)
    result.add_argument("--output-dir", type=Path, default=Path("."), metavar="DIR")
    result.add_argument(
        "--combined",
        type=Path,
        metavar="BAM",
        help="write all selected genes to one BAM",
    )
    result.add_argument(
        "--index", action="store_true", help="create BAI indexes for output BAMs"
    )
    result.add_argument("--force", action="store_true", help="replace existing outputs")
    result.add_argument(
        "--dry-run", action="store_true", help="resolve regions without extracting BAMs"
    )
    result.add_argument(
        "--manifest", type=Path, metavar="TSV", help="write a run manifest"
    )
    result.add_argument(
        "--summary", type=Path, metavar="TSV", help="write per-gene coverage QC"
    )
    result.add_argument(
        "--write-bed", type=Path, metavar="BED", help="write resolved padded regions"
    )
    result.add_argument(
        "--plot",
        type=Path,
        metavar="SVG",
        help="write a binned coverage plot (exactly one gene)",
    )
    result.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    return result


def requested_tokens(args: argparse.Namespace) -> List[str]:
    tokens = list(args.gene) + list(args.genes)
    for path in args.gene_file:
        if not path.is_file():
            raise NanoFetchError(f"Gene file not found: {path}")
        for line in path.read_text(encoding="utf-8").splitlines():
            content = line.partition("#")[0].replace(",", " ")
            tokens.extend(content.split())
    for panel in args.panel:
        tokens.extend(load_panel(panel))
    return tokens


def print_panels() -> None:
    descriptions = panel_descriptions()
    for name in available_panels():
        print(f"{name}\t{descriptions.get(name, '')}")


def write_manifest(
    path: Path,
    assembly: str,
    source: Path,
    results,
    force: bool = False,
) -> None:
    if path.exists() and not force:
        raise NanoFetchError(
            f"Manifest already exists: {path}. Use --force to replace it."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = annotation_metadata()[assembly]
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "input",
                    "assembly",
                    "annotation",
                    "gene",
                    "regions",
                    "alignments",
                    "output",
                    "index",
                ]
            )
            for result in results:
                regions = ",".join(
                    f"{region.contig}:{region.start + 1}-{region.end}"
                    for region in result.regions
                )
                writer.writerow(
                    [
                        source,
                        ASSEMBLY_LABELS[assembly],
                        metadata["label"],
                        result.symbol,
                        regions,
                        result.alignments,
                        result.output,
                        result.index or "",
                    ]
                )
        os.replace(str(temporary), str(path))
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _input_paths(path: Path) -> Tuple[Path, ...]:
    """Return an alignment and conventional paths for its indexes."""
    resolved = path.resolve()
    candidates = (path, resolved)
    for bam in (path, resolved):
        candidates += (
            Path(str(bam) + ".bai"),
            Path(str(bam) + ".csi"),
            Path(str(bam) + ".crai"),
            bam.with_suffix(".bai"),
            bam.with_suffix(".csi"),
            bam.with_suffix(".crai"),
        )
    return tuple(dict.fromkeys(candidate.resolve() for candidate in candidates))


def _protected_paths(args: argparse.Namespace) -> set:
    protected = set(_input_paths(args.input))
    protected.update(path.resolve() for path in args.gene_file)
    if args.reference:
        protected.update(
            {
                args.reference.resolve(),
                Path(str(args.reference) + ".fai").resolve(),
                Path(str(args.reference) + ".gzi").resolve(),
            }
        )
    return protected


def _planned_mutations(args: argparse.Namespace, symbols: Sequence[str]):
    """Return every path the run may create, replace, or remove."""
    planned = []
    bams = []
    if not args.dry_run:
        bams = (
            [args.combined]
            if args.combined
            else [args.output_dir / output_name(symbol) for symbol in symbols]
        )
    for index, bam in enumerate(bams):
        label = "combined output" if args.combined else f"BAM for {symbols[index]}"
        planned.append((label, bam))
        index_label = "combined output" if args.combined else symbols[index]
        planned.append((f"BAI for {index_label}", Path(str(bam) + ".bai")))
        planned.append((f"CSI for {index_label}", Path(str(bam) + ".csi")))
    if args.plot:
        planned.append(("coverage plot", args.plot))
    if args.manifest:
        planned.append(("manifest", args.manifest))
    if args.summary:
        planned.append(("summary", args.summary))
    if args.write_bed:
        planned.append(("BED", args.write_bed))
    return planned


def validate_write_plan(args: argparse.Namespace, symbols: Sequence[str]) -> None:
    """Reject destructive or ambiguous path combinations before writing."""
    protected = _protected_paths(args)
    seen = {}
    planned = _planned_mutations(args, symbols)
    for label, path in planned:
        resolved = path.resolve()
        if resolved in protected:
            raise NanoFetchError(
                f"Refusing to replace an input or its index as the {label}: {path}"
            )
        if resolved in seen:
            previous = seen[resolved]
            raise NanoFetchError(
                f"Output path conflict: {previous} and {label} both use {path}."
            )
        seen[resolved] = label
        if path.is_dir():
            raise NanoFetchError(f"Output path is a directory: {path}")

    if not args.force:
        for label, path in planned:
            if path.exists():
                raise NanoFetchError(
                    f"{label.capitalize()} already exists: {path}. "
                    "Use --force to replace it."
                )


def run(argv: Sequence[str] = None) -> int:
    args = parser().parse_args(argv)
    if args.list_panels:
        print_panels()
        return 0
    if args.input is None:
        raise NanoFetchError("An input BAM or CRAM is required.")
    if not args.input.is_file():
        raise NanoFetchError(f"Input does not exist: {args.input}")
    if args.reference and not args.reference.is_file():
        raise NanoFetchError(f"Reference FASTA not found: {args.reference}")
    if args.combined and args.combined.suffix.lower() != ".bam":
        raise NanoFetchError("--combined output must use the .bam extension.")
    if args.summary and args.summary.suffix.lower() != ".tsv":
        raise NanoFetchError("--summary output must use the .tsv extension.")
    if args.write_bed and args.write_bed.suffix.lower() != ".bed":
        raise NanoFetchError("--write-bed output must use the .bed extension.")
    if args.manifest and args.dry_run:
        raise NanoFetchError("--manifest cannot be combined with --dry-run.")
    if args.summary and args.dry_run:
        raise NanoFetchError("--summary cannot be combined with --dry-run.")

    tokens = requested_tokens(args)
    mode = "rc" if args.input.suffix.lower() == ".cram" else "rb"
    open_options = {"threads": args.threads}
    if mode == "rc" and args.reference:
        open_options["reference_filename"] = str(args.reference)

    with pysam.AlignmentFile(str(args.input), mode, **open_options) as source:
        validate_input(source, args.input)
        references = list(zip(source.references, source.lengths))
        assembly = choose_assembly(args.genome, references)
        if not tokens:
            raise NanoFetchError("Select at least one gene or use --panel NAME.")
        symbols = unique_symbols(tokens, assembly)
        if args.plot and len(symbols) != 1:
            raise NanoFetchError("--plot requires exactly one resolved gene.")
        if args.plot and args.dry_run:
            raise NanoFetchError("--plot cannot be combined with --dry-run.")
        if args.plot and args.plot.suffix.lower() != ".svg":
            raise NanoFetchError("--plot output must use the .svg extension.")
        intervals_by_gene = {
            symbol: resolve_gene(symbol, assembly) for symbol in symbols
        }
        header_contigs = dict(references)
        regions_by_gene = {
            symbol: padded_regions(intervals, header_contigs, args.padding)
            for symbol, intervals in intervals_by_gene.items()
        }
        validate_write_plan(args, symbols)
        metadata = annotation_metadata()[assembly]
        print(f"Genome: {ASSEMBLY_LABELS[assembly]}", file=sys.stderr)
        print(f"Annotation: {metadata['label']}", file=sys.stderr)

        if args.dry_run:
            for symbol in symbols:
                regions = regions_by_gene[symbol]
                region_text = ", ".join(
                    f"{region.contig}:{region.start + 1}-{region.end}"
                    for region in regions
                )
                output = args.combined or args.output_dir / output_name(symbol)
                print(f"{symbol}\t{region_text}\t{output}")
            if args.write_bed:
                write_bed(args.write_bed, regions_by_gene, force=args.force)
                print(f"Wrote {args.write_bed}")
            return 0

        results = []
        if args.combined:
            regions = combined_regions(regions_by_gene.values(), header_contigs)
            result = extract_regions(
                source=source,
                symbol=",".join(symbols),
                regions=regions,
                output=args.combined,
                include_supplementary=args.include_supplementary,
                include_secondary=args.include_secondary,
                threads=args.threads,
                index_output=args.index,
                force=args.force,
            )
            results.append(result)
            index_note = f" + {result.index.name}" if result.index else ""
            print(f"Wrote {result.output} ({result.alignments} alignments){index_note}")
        else:
            for symbol in symbols:
                result = extract_gene(
                    source=source,
                    symbol=symbol,
                    intervals=intervals_by_gene[symbol],
                    output_dir=args.output_dir,
                    padding=args.padding,
                    include_supplementary=args.include_supplementary,
                    include_secondary=args.include_secondary,
                    threads=args.threads,
                    index_output=args.index,
                    force=args.force,
                )
                results.append(result)
                index_note = f" + {result.index.name}" if result.index else ""
                message = f"Wrote {result.output} ({result.alignments} alignments)"
                print(message + index_note)

        if args.plot:
            coverage_summary = summarize_coverage(
                bam=results[0].output,
                symbol=symbols[0],
                intervals=intervals_by_gene[symbols[0]],
                regions=results[0].regions,
                assembly=ASSEMBLY_LABELS[assembly],
                annotation=metadata["label"],
                source=args.input.name,
            )
            write_coverage_svg(args.plot, coverage_summary, force=args.force)
            print(f"Wrote {args.plot}")

        if args.write_bed:
            write_bed(args.write_bed, regions_by_gene, force=args.force)
            print(f"Wrote {args.write_bed}")

        if args.summary:
            qc_summaries = [
                summarize_gene(
                    source,
                    symbol,
                    intervals_by_gene[symbol],
                    header_contigs,
                    include_supplementary=args.include_supplementary,
                    include_secondary=args.include_secondary,
                )
                for symbol in symbols
            ]
            write_summary(args.summary, qc_summaries, force=args.force)
            print(f"Wrote {args.summary}")

    if args.manifest:
        write_manifest(args.manifest, assembly, args.input, results, force=args.force)
        print(f"Wrote {args.manifest}")
    return 0


def main() -> None:
    try:
        raise SystemExit(run())
    except NanoFetchError as error:
        print(f"nanofetch: error: {error}", file=sys.stderr)
        raise SystemExit(2)
    except (OSError, ValueError) as error:
        print(f"nanofetch: error: {error}", file=sys.stderr)
        raise SystemExit(2)
