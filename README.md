# NanoFetch

Ridiculously easy gene-region BAM extraction.

```console
nanofetch tumor.bam EGFR --index
Genome: GRCh38
Annotation: GENCODE 50 (GRCh38.p14)
Wrote EGFR.bam (18432 alignments) + EGFR.bam.bai
```

`nanofetch` turns gene symbols into reproducible genomic intervals, adds padding,
and extracts overlapping alignments from a coordinate-sorted BAM. It supports
GRCh37, GRCh38, and the telomere-to-telomere `T2T-CHM13v2.0` (`hs1`) assembly.

## Install

Python 3.9 or newer is required.

```bash
python -m pip install nanofetch
```

For an isolated command-line installation, use `pipx install nanofetch`.

For local development:

```bash
git clone https://github.com/AlexanderM-M/Nanofetch.git
cd Nanofetch
python -m pip install -e '.[test]'
pytest
```

`pysam` is the only runtime dependency. A separate `samtools` installation is not
required.

## Usage

One gene, using the default 1 Mb padding on both sides:

```bash
nanofetch tumor.bam EGFR
```

Several genes:

```bash
nanofetch tumor.bam --genes EGFR CDK4 PDGFRA CDKN2A MET
```

Compact families are accepted:

```bash
nanofetch tumor.bam CDKN2A/B NTRK1/2/3
```

Use the built-in CNS convenience panel:

```bash
nanofetch tumor.bam --panel cns --index --output-dir cns-regions
```

Use T2T-CHM13 explicitly, or let `nanofetch` detect it from chromosome lengths:

```bash
nanofetch t2t-aligned.bam EGFR --genome t2t
nanofetch t2t-aligned.bam --panel cns --genome hs1
```

Inspect the resolved regions without writing BAMs:

```console
nanofetch tumor.bam EGFR --dry-run
Genome: GRCh38
Annotation: GENCODE 50 (GRCh38.p14)
EGFR   chr7:54018820-56211628   EGFR.bam
```

Generate a reproducibility manifest:

```bash
nanofetch tumor.bam EGFR MET --index --manifest regions.tsv
```

### SVG coverage plot

Create a self-contained coverage overview for one extracted gene:

```bash
nanofetch tumor.bam EGFR --plot EGFR.coverage.svg
```

The plot shows mean aligned-base depth in fixed-width bins, the annotated gene
body, requested extraction span, assembly and annotation provenance, and counts
for primary, supplementary, and secondary alignments. It visualizes the output
BAM, so `--include-supplementary` and `--include-secondary` are reflected in the
plot. SVG generation has no dependency beyond NanoFetch itself.

![Example NanoFetch coverage plot](https://raw.githubusercontent.com/AlexanderM-M/Nanofetch/main/docs/example-coverage.svg)

Run `nanofetch --help` for all options.

## What is selected?

An alignment is selected when its aligned reference span overlaps the padded gene
interval. It does not need to span the entire region.

By default, secondary (`0x100`) and supplementary (`0x800`) records are excluded.
Use `--include-secondary` and `--include-supplementary` to retain records of those
types that overlap the region.

`--include-supplementary` does **not** retrieve distant segments merely because a
different segment from the same read overlaps the requested gene. That operation
requires a read-name-based second pass and is intentionally outside the v0.1
semantics.

Overlapping padded intervals for multiple copies of the same symbol are merged, so
an alignment is not written twice. Separate gene output BAMs may contain the same
alignment when their padded regions overlap; this is expected.

## Genome detection and contig names

Automatic detection compares exact primary-chromosome lengths in the BAM header.
It will stop rather than guess when the build is ambiguous. The explicit spellings
below are accepted:

| Assembly | Accepted names | Bundled annotation |
|---|---|---|
| GRCh37 | `grch37`, `hg19` | GENCODE 50lift37 |
| GRCh38 | `grch38`, `hg38` | GENCODE 50 |
| T2T-CHM13v2.0 | `t2t`, `hs1`, `chm13`, `chm13v2.0` | RefSeq/Liftoff v5.3 |

Both `chr7` and `7` naming styles are supported. T2T GenBank names such as
`CP068271.2` and RefSeq names such as `NC_060931.1` are also mapped to their
chromosomes.

## Input and output safety

The input must be a coordinate-sorted BAM with a BAI or CSI index. `nanofetch`
does not modify the input or create an input index automatically.

Existing output BAMs are protected unless `--force` is supplied. Outputs are
written to a temporary file and moved into place only after writing succeeds.
When an output is replaced, any stale BAI or CSI is removed. `--index` creates a
fresh BAI.

## Built-in CNS panel

The `cns` panel currently contains:

```text
EGFR PDGFRA CDK4 MDM2 MDM4 MET CDKN2A CDKN2B NF1 PTEN TERT
BRAF FGFR3 NTRK1 NTRK2 NTRK3 ALK
```

This is a transparent convenience set, not a validated clinical assay. Use
`nanofetch --list-panels` to inspect installed panels.

## Annotation provenance

The package works offline. It bundles compact gene-level tables generated from
pinned upstream files; it does not bundle complete GTF/GFF files.

- GRCh38: [GENCODE release 50](https://www.gencodegenes.org/human/)
- GRCh37: [GENCODE 50lift37](https://www.gencodegenes.org/human/grch37_mapped_releases.html)
- T2T: [T2T-CHM13v2.0 RefSeq/Liftoff v5.3](https://github.com/marbl/CHM13#gene-annotation)

Source URLs, SHA-256 digests, labels, and record counts are stored in
`src/nanofetch/data/annotations.json`. The deterministic build script is
`scripts/build_annotations.py`.

## Scope

Version 0.1 intentionally writes one BAM per gene. Potential later additions
include custom BED/GTF annotations, CRAM output, mate retrieval, all-segment
retrieval, and combined multi-gene BAMs.

## License

MIT. T2T consortium data are released under CC0; see the upstream
[CHM13 repository](https://github.com/marbl/CHM13#data-reuse-and-license).
