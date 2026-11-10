# Troubleshooting

[Back to the README](../README.md) · [Usage guide](usage.md)

## Command not found or module not found

Install NanoFetch in the Python environment you intend to use, then check it:

```bash
python -m pip install nanofetch
python -m nanofetch --version
```

If `python -m nanofetch` works but `nanofetch` does not, activate your virtual
environment or add that environment's scripts directory to your `PATH`.
For a repository checkout, follow the editable installation in
[Contributing](../CONTRIBUTING.md).

## Input is not indexed or not coordinate-sorted

NanoFetch needs coordinate-sorted input with a matching index. It does not sort
or index the input automatically. BAM accepts BAI or CSI; CRAM uses CRAI.

For an already coordinate-sorted BAM, create its index using the installed
`pysam` dependency:

```bash
python -c 'import pysam; pysam.index("tumor.bam")'
```

For unsorted BAM input, write a separate sorted copy and index that copy:

```bash
python -c 'import pysam; pysam.sort("-o", "tumor.sorted.bam", "tumor.bam"); pysam.index("tumor.sorted.bam")'
nanofetch tumor.sorted.bam EGFR --index
```

Choose a new sorted-output filename if `tumor.sorted.bam` already exists.
The `--index` option on NanoFetch creates **output** indexes; it does not repair
a missing input index.

## Genome assembly cannot be detected

Detection compares exact primary-chromosome lengths in the alignment header.
A subset of chromosomes or nonstandard contig names may leave it ambiguous.
If you know the assembly used for alignment, supply it explicitly:

```bash
nanofetch tumor.bam EGFR --genome grch38 --dry-run
```

Use the actual alignment reference: `grch37`/`hg19`, `grch38`/`hg38`, or
`t2t`/`hs1`. Selecting a build does not convert coordinates or realign reads.
If the header conflicts with the requested build, verify the original reference.
See [assembly and contig support](usage.md#genome-detection-and-contig-names).

## Gene not found or alias is ambiguous

Check the spelling and the selected assembly. Symbols are case-insensitive;
NanoFetch suggests close matches when available. For an ambiguous alias, use
one of the canonical symbols listed in the error. A symbol may be absent from
the bundled annotation for a particular assembly; inspect the
[annotation provenance](annotations.md) when comparing annotation releases.

Gene files accept whitespace or commas between symbols and `#` comments.
Preview a panel with:

```bash
nanofetch tumor.bam --gene-file my-panel.txt --dry-run
```

## CRAM cannot find its reference

Supply the FASTA used to create the CRAM, with its FASTA index alongside it:

```bash
nanofetch tumor.cram EGFR --reference GRCh38.fa --index
```

Use a matching reference, not just a similarly named assembly file. CRAM also
needs its own CRAI index. Offline annotations do not replace the reference
sequence needed to decode CRAM.

## Output contains zero alignments

First inspect the resolved region with `--dry-run`. Check the assembly, selected
gene, and whether the input has alignments there. Default padding is 1 Mb on
each side; `--padding 0` limits extraction to the gene interval.

Secondary and supplementary alignments are excluded by default. To retain
those records when they overlap the region:

```bash
nanofetch tumor.bam EGFR --include-secondary --include-supplementary \
  --output-dir results/all-alignments --index
```

These flags do not fetch distant segments or mates. A zero-alignment BAM can
be a valid result; it does not itself imply a gene is biologically absent.

## Output already exists or appears in an unexpected directory

Use a separate `--output-dir` for each sample. It controls per-gene BAMs and
their indexes only. Combined BAMs and reports use the exact paths you pass;
see [output locations](usage.md#where-outputs-go) for examples.

Use `--force` only to replace outputs intentionally. It cannot override input
protection or resolve two outputs targeting the same path. A successful dry run
does not guarantee planned BAM filenames are available for extraction.

## Still stuck?

[Open a bug report](https://github.com/AlexanderM-M/Nanofetch/issues/new?template=bug_report.yml)
with `nanofetch --version`, your Python version and operating system, the command,
the assembly, and the full error text. A small synthetic reproducer is helpful;
do not upload patient data or private alignment files.
