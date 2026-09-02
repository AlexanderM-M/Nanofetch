from pathlib import Path

import pysam
import pytest
from conftest import _alignment, make_bam

from nanofetch.annotations import resolve_gene
from nanofetch.cli import run
from nanofetch.errors import NanoFetchError
from nanofetch.extract import (
    extract_gene,
    iter_region_alignments,
    padded_regions,
    validate_input,
)
from nanofetch.models import GeneInterval, Region


def test_padding_uses_zero_based_half_open_coordinates():
    interval = resolve_gene("EGFR", "grch38")
    regions = padded_regions(interval, {"chr7": 159345973}, 1_000_000)
    assert [(r.contig, r.start, r.end) for r in regions] == [
        ("chr7", 54018819, 56211628)
    ]


def test_extract_primary_and_create_index(grch38_bam, tmp_path):
    output_dir = tmp_path / "result"
    with pysam.AlignmentFile(grch38_bam, "rb") as source:
        validate_input(source, grch38_bam)
        result = extract_gene(
            source, "EGFR", resolve_gene("EGFR", "grch38"), output_dir,
            padding=0, index_output=True,
        )
    assert result.alignments == 1
    assert result.output.exists()
    assert result.index.exists()
    with pysam.AlignmentFile(result.output, "rb") as output:
        assert [record.query_name for record in output] == ["primary"]
        assert output.header.to_dict()["PG"][-1]["PN"] == "nanofetch"


def test_supplementary_and_secondary_are_opt_in(grch38_bam, tmp_path):
    with pysam.AlignmentFile(grch38_bam, "rb") as source:
        result = extract_gene(
            source, "EGFR", resolve_gene("EGFR", "grch38"), tmp_path / "all",
            padding=0, include_supplementary=True, include_secondary=True,
        )
    assert result.alignments == 3


def test_existing_output_is_protected(grch38_bam, tmp_path):
    output_dir = tmp_path / "result"
    output_dir.mkdir()
    (output_dir / "EGFR.bam").write_bytes(b"do not replace")
    with pysam.AlignmentFile(grch38_bam, "rb") as source:
        with pytest.raises(NanoFetchError, match="already exists"):
            extract_gene(
                source, "EGFR", resolve_gene("EGFR", "grch38"), output_dir,
            )
    assert (output_dir / "EGFR.bam").read_bytes() == b"do not replace"


def test_cli_auto_detects_and_writes_manifest(grch38_bam, tmp_path, capsys):
    output_dir = tmp_path / "cli"
    manifest = tmp_path / "manifest.tsv"
    assert run([
        str(grch38_bam), "EGFR", "--padding", "0", "--output-dir",
        str(output_dir), "--index", "--manifest", str(manifest),
    ]) == 0
    assert (output_dir / "EGFR.bam").exists()
    assert (output_dir / "EGFR.bam.bai").exists()
    assert "T2T" not in capsys.readouterr().err
    assert "GRCh38" in manifest.read_text()


def test_t2t_accession_named_bam(tmp_path):
    input_bam = make_bam(tmp_path / "t2t.bam", assembly="t2t", accession=True)
    output_dir = tmp_path / "t2t-output"
    assert run([
        str(input_bam), "EGFR", "--padding", "0", "--output-dir", str(output_dir)
    ]) == 0
    with pysam.AlignmentFile(output_dir / "EGFR.bam", "rb") as output:
        assert [record.query_name for record in output] == ["primary"]


def test_unindexed_input_has_actionable_error(tmp_path):
    bam = make_bam(tmp_path / "unindexed.bam")
    Path(str(bam) + ".bai").unlink()
    with pysam.AlignmentFile(bam, "rb") as source:
        with pytest.raises(NanoFetchError, match="samtools index"):
            validate_input(source, bam)


def test_force_cannot_replace_the_input_bam(tmp_path):
    input_bam = make_bam(tmp_path / "EGFR.bam")
    original = input_bam.read_bytes()

    with pytest.raises(NanoFetchError, match="input or its index"):
        run([
            str(input_bam), "EGFR", "--padding", "0",
            "--output-dir", str(tmp_path), "--force",
        ])

    assert input_bam.read_bytes() == original
    with pysam.AlignmentFile(input_bam, "rb") as source:
        source.check_index()


@pytest.mark.parametrize("target_kind", ["bam", "index"])
def test_manifest_cannot_replace_input_or_its_index(tmp_path, target_kind):
    input_bam = make_bam(tmp_path / "input.bam")
    input_index = Path(str(input_bam) + ".bai")
    target = input_bam if target_kind == "bam" else input_index
    original = target.read_bytes()

    with pytest.raises(NanoFetchError, match="input or its index"):
        run([
            str(input_bam), "EGFR", "--padding", "0",
            "--output-dir", str(tmp_path / "output"),
            "--manifest", str(target), "--force",
        ])

    assert target.read_bytes() == original
    assert not (tmp_path / "output" / "EGFR.bam").exists()


def test_all_output_collisions_are_checked_before_extraction(grch38_bam, tmp_path):
    output_dir = tmp_path / "output"
    shared = tmp_path / "report.svg"

    with pytest.raises(NanoFetchError, match="Output path conflict"):
        run([
            str(grch38_bam), "EGFR", "--padding", "0",
            "--output-dir", str(output_dir),
            "--plot", str(shared), "--manifest", str(shared),
        ])

    assert not (output_dir / "EGFR.bam").exists()
    assert not shared.exists()


def test_manifest_is_protected_and_force_can_replace_it(grch38_bam, tmp_path):
    output_dir = tmp_path / "output"
    manifest = tmp_path / "regions.tsv"
    manifest.write_text("keep me", encoding="utf-8")

    with pytest.raises(NanoFetchError, match="Manifest already exists"):
        run([
            str(grch38_bam), "EGFR", "--padding", "0",
            "--output-dir", str(output_dir), "--manifest", str(manifest),
        ])
    assert manifest.read_text(encoding="utf-8") == "keep me"
    assert not (output_dir / "EGFR.bam").exists()

    assert run([
        str(grch38_bam), "EGFR", "--padding", "0",
        "--output-dir", str(output_dir), "--manifest", str(manifest),
        "--force",
    ]) == 0
    assert manifest.read_text(encoding="utf-8").startswith(
        "input\tassembly\tannotation\tgene\t"
    )


def test_dry_run_reports_the_sanitized_output_name(tmp_path, capsys):
    input_bam = make_bam(
        tmp_path / "t2t.bam", assembly="t2t", contig="chr5"
    )

    assert run([
        str(input_bam), "PCDHA@", "--genome", "t2t",
        "--padding", "0", "--dry-run",
        "--output-dir", str(tmp_path / "output"),
    ]) == 0

    output = capsys.readouterr().out
    assert "PCDHA_.bam" in output
    assert "PCDHA@.bam" not in output


def test_sanitized_collision_is_found_before_any_gene_is_written(tmp_path):
    input_bam = make_bam(
        tmp_path / "t2t.bam", assembly="t2t", contig="chr5"
    )
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    existing = output_dir / "PCDHA_.bam"
    existing.write_bytes(b"keep me")

    with pytest.raises(NanoFetchError, match="already exists"):
        run([
            str(input_bam), "PCDHB@", "PCDHA@", "--genome", "t2t",
            "--padding", "0", "--output-dir", str(output_dir),
        ])

    assert existing.read_bytes() == b"keep me"
    assert not (output_dir / "PCDHB_.bam").exists()


def test_combined_output_contains_multiple_genes_once(tmp_path):
    intervals = {
        symbol: resolve_gene(symbol, "grch38")[0] for symbol in ("EGFR", "MET")
    }
    input_bam = tmp_path / "input.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr7", "LN": 159345973}],
    }
    records = sorted(
        (
            _alignment(symbol.lower(), interval.start - 1)
            for symbol, interval in intervals.items()
        ),
        key=lambda record: record.reference_start,
    )
    with pysam.AlignmentFile(input_bam, "wb", header=header) as output:
        for record in records:
            output.write(record)
    pysam.index(str(input_bam))
    combined = tmp_path / "panel.bam"

    assert run([
        str(input_bam), "EGFR", "MET", "--padding", "0",
        "--combined", str(combined), "--index",
    ]) == 0

    with pysam.AlignmentFile(combined, "rb") as output:
        assert [record.query_name for record in output] == [
            record.query_name for record in records
        ]
    assert Path(str(combined) + ".bai").is_file()
    assert not (tmp_path / "EGFR.bam").exists()


def test_region_union_deduplicates_a_record_spanning_two_regions(tmp_path):
    input_bam = tmp_path / "spanning.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(input_bam, "wb", header=header) as output:
        output.write(_alignment("spanning", 50))
    pysam.index(str(input_bam))

    with pysam.AlignmentFile(input_bam, "rb") as source:
        records = list(
            iter_region_alignments(
                source, (Region("chr1", 0, 100), Region("chr1", 120, 160))
            )
        )
    assert [record.query_name for record in records] == ["spanning"]


def test_gene_file_and_bed_only_dry_run(grch38_bam, tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("# custom panel\nEGFR, EGFR  # duplicate\n", encoding="utf-8")
    bed = tmp_path / "regions.bed"

    assert run([
        str(grch38_bam), "--gene-file", str(genes), "--padding", "0",
        "--write-bed", str(bed), "--dry-run",
    ]) == 0

    interval = resolve_gene("EGFR", "grch38")[0]
    assert bed.read_text(encoding="utf-8") == (
        f"chr7\t{interval.start - 1}\t{interval.end}\tEGFR\n"
    )
    assert not (tmp_path / "EGFR.bam").exists()


def test_summary_reports_gene_body_coverage(grch38_bam, tmp_path):
    summary = tmp_path / "coverage.tsv"
    assert run([
        str(grch38_bam), "EGFR", "--padding", "0",
        "--output-dir", str(tmp_path / "output"), "--summary", str(summary),
    ]) == 0

    header, row = [line.split("\t") for line in summary.read_text().splitlines()]
    values = dict(zip(header, row))
    assert values["gene"] == "EGFR"
    assert values["alignments"] == "1"
    assert values["mean_mapq"] == "60.00"
    assert float(values["mean_depth"]) > 0
    assert float(values["covered_1x_pct"]) > 0
    assert values["covered_10x_pct"] == "0.000"


def test_cram_input_is_supported(grch38_bam, tmp_path):
    cram = tmp_path / "input.cram"
    with pysam.AlignmentFile(grch38_bam, "rb") as source:
        with pysam.AlignmentFile(
            cram, "wc", header=source.header, format_options=[b"no_ref=1"]
        ) as output:
            for alignment in source:
                output.write(alignment)
    pysam.index(str(cram))

    output_dir = tmp_path / "cram-output"
    assert run([
        str(cram), "EGFR", "--padding", "0", "--output-dir", str(output_dir)
    ]) == 0
    with pysam.AlignmentFile(output_dir / "EGFR.bam", "rb") as output:
        assert [record.query_name for record in output] == ["primary"]


def test_cram_reference_is_used(tmp_path, monkeypatch):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr7\n" + "A" * 1000 + "\n", encoding="ascii")
    pysam.faidx(str(reference))
    cram = tmp_path / "input.cram"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr7", "LN": 1000}],
    }
    with pysam.AlignmentFile(
        cram, "wc", header=header, reference_filename=str(reference)
    ) as output:
        output.write(_alignment("primary", 100))
    pysam.index(str(cram))
    interval = GeneInterval("EGFR", "gene", "chr7", 101, 200, "protein_coding")
    monkeypatch.setattr("nanofetch.cli.choose_assembly", lambda *_: "grch38")
    monkeypatch.setattr("nanofetch.cli.resolve_gene", lambda *_: (interval,))

    output_dir = tmp_path / "output"
    assert run(
        [
            str(cram),
            "EGFR",
            "--reference",
            str(reference),
            "--padding",
            "0",
            "--output-dir",
            str(output_dir),
        ]
    ) == 0
    with pysam.AlignmentFile(output_dir / "EGFR.bam", "rb") as output:
        assert [record.query_name for record in output] == ["primary"]
