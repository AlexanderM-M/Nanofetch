from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


@dataclass(frozen=True)
class GeneInterval:
    symbol: str
    gene_id: str
    contig: str
    start: int
    end: int
    biotype: str
    aliases: Tuple[str, ...] = ()


@dataclass(frozen=True)
class Region:
    contig: str
    start: int  # zero-based, inclusive
    end: int  # zero-based, exclusive


@dataclass(frozen=True)
class ExtractionResult:
    symbol: str
    output: Path
    regions: Tuple[Region, ...]
    alignments: int
    index: Optional[Path] = None
