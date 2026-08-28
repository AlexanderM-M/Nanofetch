from pathlib import Path

import pysam
import pytest

from nanofetch.annotations import resolve_gene
from nanofetch.cli import run
from nanofetch.errors import NanoFetchError
from nanofetch.extract import extract_gene, padded_regions, validate_input

from conftest import make_bam


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

    with pytest.raises(NanoFetchError, match="input BAM or its index"):
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

    with pytest.raises(NanoFetchError, match="input BAM or its index"):
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
