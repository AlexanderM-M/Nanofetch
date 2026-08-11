import pytest

from bamregions.annotations import annotation_metadata, resolve_gene, unique_symbols
from bamregions.errors import GeneNotFoundError
from bamregions.panels import load_panel


@pytest.mark.parametrize("assembly", ["grch37", "grch38", "t2t"])
def test_cns_panel_resolves_on_every_assembly(assembly):
    for symbol in load_panel("cns"):
        assert resolve_gene(symbol, assembly)


def test_t2t_egfr_coordinates_and_provenance():
    egfr = resolve_gene("egfr", "t2t")
    assert len(egfr) == 1
    assert (egfr[0].contig, egfr[0].start, egfr[0].end) == (
        "chr7", 55178937, 55372056
    )
    assert "v5.3" in annotation_metadata()["t2t"]["label"]


def test_compact_gene_notation():
    assert unique_symbols(["CDKN2A/B", "NTRK1/2/3"], "grch38") == [
        "CDKN2A", "CDKN2B", "NTRK1", "NTRK2", "NTRK3"
    ]


def test_unknown_gene_has_clean_error():
    with pytest.raises(GeneNotFoundError, match="not present"):
        resolve_gene("DEFINITELY_NOT_A_GENE", "grch38")

