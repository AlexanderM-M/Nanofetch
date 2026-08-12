import csv
import difflib
import gzip
import json
from functools import lru_cache
from importlib import resources
from typing import Dict, List, Sequence, Tuple

from .errors import NanoFetchError, GeneNotFoundError
from .models import GeneInterval


DATA_PACKAGE = "nanofetch.data"


@lru_cache(maxsize=None)
def annotation_metadata() -> Dict[str, dict]:
    resource = resources.files(DATA_PACKAGE).joinpath("annotations.json")
    with resource.open("r", encoding="utf-8") as handle:
        return json.load(handle)


@lru_cache(maxsize=3)
def load_annotation(assembly: str) -> Dict[str, Tuple[GeneInterval, ...]]:
    resource = resources.files(DATA_PACKAGE).joinpath(f"{assembly}.genes.tsv.gz")
    genes: Dict[str, List[GeneInterval]] = {}
    with resource.open("rb") as raw, gzip.open(raw, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            aliases = tuple(filter(None, row["aliases"].split(",")))
            interval = GeneInterval(
                symbol=row["symbol"], gene_id=row["gene_id"],
                contig=row["contig"], start=int(row["start"]),
                end=int(row["end"]), biotype=row["biotype"], aliases=aliases,
            )
            keys = {row["symbol"].upper(), *(alias.upper() for alias in aliases)}
            for key in keys:
                genes.setdefault(key, []).append(interval)
    return {key: tuple(value) for key, value in genes.items()}


def resolve_gene(symbol: str, assembly: str) -> Tuple[GeneInterval, ...]:
    genes = load_annotation(assembly)
    key = symbol.upper()
    if key in genes:
        matches = genes[key]
        exact = tuple(interval for interval in matches if interval.symbol.upper() == key)
        if exact:
            return exact
        canonical_symbols = sorted({interval.symbol for interval in matches})
        if len(canonical_symbols) > 1:
            raise NanoFetchError(
                f"Gene alias {symbol!r} is ambiguous in the {assembly} annotation: "
                f"{', '.join(canonical_symbols)}. Use a canonical symbol."
            )
        return matches
    suggestions = difflib.get_close_matches(key, genes.keys(), n=3, cutoff=0.72)
    hint = f" Did you mean {', '.join(suggestions)}?" if suggestions else ""
    raise GeneNotFoundError(
        f"Gene {symbol!r} is not present in the {assembly} annotation.{hint}"
    )


def expand_gene_token(token: str, assembly: str) -> List[str]:
    """Expand compact forms such as CDKN2A/B and NTRK1/2/3."""
    if "/" not in token:
        return [token]
    parts = token.split("/")
    if not all(parts):
        return [token]
    genes = load_annotation(assembly)
    first = parts[0]
    expanded = [first]
    tail_is_digit = all(part.isdigit() for part in parts[1:])
    tail_is_alpha = all(part.isalpha() for part in parts[1:])
    if tail_is_digit:
        base = first.rstrip("0123456789")
    elif tail_is_alpha:
        base = first.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
    else:
        base = ""
    for part in parts[1:]:
        candidate = part if part.upper() in genes else base + part
        expanded.append(candidate)
    return expanded


def unique_symbols(tokens: Sequence[str], assembly: str) -> List[str]:
    result: List[str] = []
    seen = set()
    for token in tokens:
        for expanded in expand_gene_token(token, assembly):
            intervals = resolve_gene(expanded, assembly)
            canonical = intervals[0].symbol
            if canonical.upper() not in seen:
                result.append(canonical)
                seen.add(canonical.upper())
    return result
