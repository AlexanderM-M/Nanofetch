# Bundled annotation data

The `*.genes.tsv.gz` files contain gene-level, 1-based inclusive intervals. They
are deterministic reductions of the upstream annotations described in
`annotations.json`.

To regenerate after downloading the exact upstream files:

```bash
python scripts/build_annotations.py \
  --grch37 gencode.v50lift37.annotation.gtf.gz \
  --grch38 gencode.v50.annotation.gtf.gz \
  --t2t chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz \
  --output src/bamregions/data
```

The CNS panel is a convenience list and is not a validated clinical assay.

