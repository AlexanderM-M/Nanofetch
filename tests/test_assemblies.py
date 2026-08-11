import pytest

from bamregions.assemblies import (
    choose_assembly,
    detect_assembly,
    normalize_assembly,
    resolve_contig,
)
from bamregions.errors import AssemblyError


@pytest.mark.parametrize(
    ("reference", "length", "expected"),
    [
        ("chr1", 249250621, "grch37"),
        ("1", 248956422, "grch38"),
        ("chr1", 248387328, "t2t"),
        ("CP068277.2", 248387328, "t2t"),
        ("NC_060925.1", 248387328, "t2t"),
    ],
)
def test_detect_assembly(reference, length, expected):
    assert detect_assembly([(reference, length)]) == expected


def test_explicit_mismatch_is_rejected():
    with pytest.raises(AssemblyError, match="matches GRCh38"):
        choose_assembly("grch37", [("chr7", 159345973)])


def test_ambiguous_build_is_rejected():
    with pytest.raises(AssemblyError, match="Could not identify"):
        detect_assembly([("custom", 123)])


@pytest.mark.parametrize("alias", ["t2t", "hs1", "chm13", "chm13v2.0", "t2t-chm13v2.0"])
def test_t2t_names(alias):
    assert normalize_assembly(alias) == "t2t"


def test_resolve_t2t_accession_contig():
    assert resolve_contig("chr7", {"CP068271.2": 160567428}) == "CP068271.2"
    assert resolve_contig("chr7", {"NC_060931.1": 160567428}) == "NC_060931.1"

