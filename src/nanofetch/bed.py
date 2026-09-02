"""BED export for resolved extraction regions."""

import os
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

from .errors import NanoFetchError
from .models import Region


def write_bed(
    path: Path, regions_by_gene: Mapping[str, Sequence[Region]], force: bool = False
) -> None:
    """Atomically write zero-based, half-open BED4 records."""
    if path.exists() and not force:
        raise NanoFetchError(f"BED already exists: {path}. Use --force to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            for symbol, regions in regions_by_gene.items():
                for region in regions:
                    handle.write(
                        f"{region.contig}\t{region.start}\t{region.end}\t{symbol}\n"
                    )
        os.replace(str(temporary), str(path))
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
