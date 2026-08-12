import argparse
import csv
import sys
from pathlib import Path
from typing import List, Sequence

import pysam

from . import __version__
from .annotations import annotation_metadata, resolve_gene, unique_symbols
from .assemblies import ASSEMBLY_LABELS, choose_assembly
from .errors import NanoFetchError
from .extract import extract_gene, padded_regions, validate_input
from .panels import available_panels, load_panel, panel_descriptions
from .plot import summarize_coverage, write_coverage_svg


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
        description="Ridiculously easy gene-region BAM extraction.",
    )
    result.add_argument("input", nargs="?", type=Path, help="coordinate-sorted, indexed BAM")
    result.add_argument("gene", nargs="*", help="one or more gene symbols")
    result.add_argument("--genes", nargs="+", default=[], metavar="GENE", help="additional gene symbols")
    result.add_argument("--panel", action="append", default=[], metavar="NAME", help="built-in panel")
    result.add_argument("--list-panels", action="store_true", help="list built-in panels and exit")
    result.add_argument(
        "--genome", default="auto", metavar="BUILD",
        help="auto, grch37/hg19, grch38/hg38, or t2t/hs1 (default: auto)",
    )
    result.add_argument("--padding", type=nonnegative_integer, default=1_000_000, metavar="BP")
    result.add_argument("--include-supplementary", action="store_true")
    result.add_argument("--include-secondary", action="store_true")
    result.add_argument("--threads", type=positive_integer, default=1)
    result.add_argument("--output-dir", type=Path, default=Path("."), metavar="DIR")
    result.add_argument("--index", action="store_true", help="create BAI indexes for output BAMs")
    result.add_argument("--force", action="store_true", help="replace existing output BAMs")
    result.add_argument("--dry-run", action="store_true", help="resolve and report without writing")
    result.add_argument("--manifest", type=Path, metavar="TSV", help="write a run manifest")
    result.add_argument(
        "--plot", type=Path, metavar="SVG",
        help="write a binned coverage plot (exactly one gene)",
    )
    result.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return result


def requested_tokens(args: argparse.Namespace) -> List[str]:
    tokens = list(args.gene) + list(args.genes)
    for panel in args.panel:
        tokens.extend(load_panel(panel))
    return tokens


def print_panels() -> None:
    descriptions = panel_descriptions()
    for name in available_panels():
        print(f"{name}\t{descriptions.get(name, '')}")


def write_manifest(path: Path, assembly: str, source: Path, results) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = annotation_metadata()[assembly]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "input", "assembly", "annotation", "gene", "regions",
            "alignments", "output", "index",
        ])
        for result in results:
            regions = ",".join(
                f"{region.contig}:{region.start + 1}-{region.end}" for region in result.regions
            )
            writer.writerow([
                source, ASSEMBLY_LABELS[assembly], metadata["label"], result.symbol,
                regions, result.alignments, result.output, result.index or "",
            ])


def run(argv: Sequence[str] = None) -> int:
    args = parser().parse_args(argv)
    if args.list_panels:
        print_panels()
        return 0
    if args.input is None:
        raise NanoFetchError("An input BAM is required.")
    if not args.input.exists():
        raise NanoFetchError(f"Input does not exist: {args.input}")

    with pysam.AlignmentFile(str(args.input), "rb", threads=args.threads) as source:
        validate_input(source, args.input)
        references = list(zip(source.references, source.lengths))
        assembly = choose_assembly(args.genome, references)
        tokens = requested_tokens(args)
        if not tokens:
            raise NanoFetchError("Select at least one gene or use --panel NAME.")
        symbols = unique_symbols(tokens, assembly)
        if args.plot and len(symbols) != 1:
            raise NanoFetchError("--plot requires exactly one resolved gene.")
        if args.plot and args.dry_run:
            raise NanoFetchError("--plot cannot be combined with --dry-run.")
        if args.plot and args.plot.suffix.lower() != ".svg":
            raise NanoFetchError("--plot output must use the .svg extension.")
        if args.plot and args.plot.exists() and not args.force:
            raise NanoFetchError(
                f"Plot already exists: {args.plot}. Use --force to replace it."
            )
        metadata = annotation_metadata()[assembly]
        print(f"Genome: {ASSEMBLY_LABELS[assembly]}", file=sys.stderr)
        print(f"Annotation: {metadata['label']}", file=sys.stderr)

        header_contigs = dict(references)
        if args.dry_run:
            for symbol in symbols:
                regions = padded_regions(resolve_gene(symbol, assembly), header_contigs, args.padding)
                region_text = ", ".join(
                    f"{region.contig}:{region.start + 1}-{region.end}" for region in regions
                )
                print(f"{symbol}\t{region_text}\t{args.output_dir / (symbol + '.bam')}")
            return 0

        # Fail before writing anything when an output collision is known.
        if not args.force:
            collisions = [args.output_dir / f"{symbol}.bam" for symbol in symbols
                          if (args.output_dir / f"{symbol}.bam").exists()]
            if collisions:
                raise NanoFetchError(
                    f"Output already exists: {collisions[0]}. Use --force to replace it."
                )

        results = []
        for symbol in symbols:
            result = extract_gene(
                source=source, symbol=symbol, intervals=resolve_gene(symbol, assembly),
                output_dir=args.output_dir, padding=args.padding,
                include_supplementary=args.include_supplementary,
                include_secondary=args.include_secondary, threads=args.threads,
                index_output=args.index, force=args.force,
            )
            results.append(result)
            index_note = f" + {result.index.name}" if result.index else ""
            print(f"Wrote {result.output} ({result.alignments} alignments){index_note}")
            if args.plot:
                summary = summarize_coverage(
                    bam=result.output,
                    symbol=symbol,
                    intervals=resolve_gene(symbol, assembly),
                    regions=result.regions,
                    assembly=ASSEMBLY_LABELS[assembly],
                    annotation=metadata["label"],
                    source=args.input.name,
                )
                write_coverage_svg(args.plot, summary, force=args.force)
                print(f"Wrote {args.plot}")

    if args.manifest:
        write_manifest(args.manifest, assembly, args.input, results)
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
