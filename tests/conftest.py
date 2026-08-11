from pathlib import Path

import pysam
import pytest

from bamregions.annotations import resolve_gene
from bamregions.assemblies import CHROMOSOME_LENGTHS


def _alignment(name, start, flag=0):
    record = pysam.AlignedSegment()
    record.query_name = name
    record.query_sequence = "A" * 100
    record.flag = flag
    record.reference_id = 0
    record.reference_start = start
    record.mapping_quality = 60
    record.cigar = ((0, 100),)
    record.query_qualities = pysam.qualitystring_to_array("I" * 100)
    return record


def make_bam(path: Path, assembly="grch38", contig="chr7", accession=False):
    interval = resolve_gene("EGFR", assembly)[0]
    chrom = interval.contig.removeprefix("chr")
    name = contig
    if accession and assembly == "t2t":
        name = "CP068271.2"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": CHROMOSOME_LENGTHS[assembly][chrom]}],
    }
    start = interval.start - 1
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        output.write(_alignment("primary", start, 0))
        output.write(_alignment("secondary", start + 1, 0x100))
        output.write(_alignment("supplementary", start + 2, 0x800))
    pysam.index(str(path))
    return path


@pytest.fixture
def grch38_bam(tmp_path):
    return make_bam(tmp_path / "input.bam")

