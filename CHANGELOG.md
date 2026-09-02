# Changelog

## 0.3.0 - unreleased

- Add deduplicated combined BAM output for multiple genes.
- Accept indexed CRAM input with an optional reference FASTA.
- Add custom gene files, BED export, and per-gene coverage QC summaries.

## 0.2.0 - 2026-09-01

- Adopt NanoFetch as the distribution, Python package, and command name.
- Add dependency-free SVG coverage plots with `--plot FILE.svg`.
- Prevent BAM, index, plot, and manifest output paths from conflicting with the
  input or with each other.
- Protect existing manifests and write them atomically.

## 0.1.0

- Extract one indexed BAM per gene with configurable padding.
- Support GRCh37, GRCh38, and T2T-CHM13v2.0.
- Bundle versioned GENCODE 50, GENCODE 50lift37, and T2T RefSeq/Liftoff v5.3
  gene intervals.
- Add assembly auto-detection, contig aliases, compact gene notation, the CNS
  panel, alignment flag controls, output indexing, dry-run, and manifests.
