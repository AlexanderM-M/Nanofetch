#!/usr/bin/env python3
"""Build compact, deterministic gene tables from pinned upstream annotations."""

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import unquote


SOURCES = {
    "grch38": {
        "label": "GENCODE 50 (GRCh38.p14)",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz",
        "format": "gtf",
    },
    "grch37": {
        "label": "GENCODE 50lift37 (GRCh37.p13)",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh37_mapping/gencode.v50lift37.annotation.gtf.gz",
        "format": "gtf",
    },
    "t2t": {
        "label": "T2T-CHM13v2.0 RefSeq/Liftoff v5.3",
        "url": "https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz",
        "format": "gff3",
    },
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gtf_attributes(text):
    result = {}
    for item in text.rstrip(";").split("; "):
        if " " in item:
            key, value = item.split(" ", 1)
            result[key] = value.strip('"')
    return result


def gff_attributes(text):
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def rows(path, file_format):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = gtf_attributes(fields[8]) if file_format == "gtf" else gff_attributes(fields[8])
            if file_format == "gtf":
                symbol = attrs.get("gene_name")
                gene_id = attrs.get("gene_id", "")
                biotype = attrs.get("gene_type", "")
                aliases = ""
            else:
                symbol = attrs.get("gene_name") or attrs.get("gene")
                gene_id = attrs.get("ID", "")
                biotype = attrs.get("gene_biotype", "")
                aliases = attrs.get("gene_synonym", "")
            if symbol:
                yield (symbol, gene_id, fields[0], int(fields[3]), int(fields[4]), biotype, aliases)


def natural_contig(value):
    name = value[3:] if value.startswith("chr") else value
    if name.isdigit():
        return (0, int(name))
    return (1, name)


def write_table(source, destination, file_format):
    records = sorted(set(rows(source, file_format)), key=lambda row: (natural_contig(row[2]), row[3], row[0]))
    with destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(["symbol", "gene_id", "contig", "start", "end", "biotype", "aliases"])
                writer.writerows(records)
    return len(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--grch37", required=True, type=Path)
    parser.add_argument("--grch38", required=True, type=Path)
    parser.add_argument("--t2t", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = {}
    for assembly, source_path in (("grch37", args.grch37), ("grch38", args.grch38), ("t2t", args.t2t)):
        info = dict(SOURCES[assembly])
        count = write_table(source_path, args.output / f"{assembly}.genes.tsv.gz", info.pop("format"))
        info.update({"source_sha256": sha256(source_path), "gene_records": count})
        metadata[assembly] = info
        print(f"{assembly}: {count} gene records")
    with (args.output / "annotations.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()

