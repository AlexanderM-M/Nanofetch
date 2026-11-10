# Annotation provenance

[Back to the README](../README.md)

The package works offline. It bundles compact gene-level tables generated from
pinned upstream files; it does not bundle complete GTF/GFF files.

- GRCh38: [GENCODE release 50](https://www.gencodegenes.org/human/)
- GRCh37: [GENCODE 50lift37](https://www.gencodegenes.org/human/grch37_mapped_releases.html)
- T2T: [T2T-CHM13v2.0 RefSeq/Liftoff v5.3](https://github.com/marbl/CHM13#gene-annotation)

Source URLs, SHA-256 digests, labels, and record counts are stored in
[`annotations.json`](../src/nanofetch/data/annotations.json). The deterministic build script is
[`scripts/build_annotations.py`](../scripts/build_annotations.py).

See the [data documentation](../src/nanofetch/data/README.md) for table formats
and rebuilding instructions, and the [usage guide](usage.md#genome-detection-and-contig-names)
for assembly aliases and contig matching.
