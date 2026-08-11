from typing import Dict, Iterable, Mapping, Sequence, Tuple

from .errors import AssemblyError


ASSEMBLY_LABELS = {
    "grch37": "GRCh37",
    "grch38": "GRCh38",
    "t2t": "T2T-CHM13v2.0",
}

# Primary chromosome lengths. Exact matches keep auto-detection conservative.
CHROMOSOME_LENGTHS: Dict[str, Dict[str, int]] = {
    "grch37": {
        "1": 249250621, "2": 243199373, "3": 198022430,
        "4": 191154276, "5": 180915260, "6": 171115067,
        "7": 159138663, "8": 146364022, "9": 141213431,
        "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392,
        "16": 90354753, "17": 81195210, "18": 78077248,
        "19": 59128983, "20": 63025520, "21": 48129895,
        "22": 51304566, "X": 155270560, "Y": 59373566,
        "M": 16571,
    },
    "grch38": {
        "1": 248956422, "2": 242193529, "3": 198295559,
        "4": 190214555, "5": 181538259, "6": 170805979,
        "7": 159345973, "8": 145138636, "9": 138394717,
        "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189,
        "16": 90338345, "17": 83257441, "18": 80373285,
        "19": 58617616, "20": 64444167, "21": 46709983,
        "22": 50818468, "X": 156040895, "Y": 57227415,
        "M": 16569,
    },
    "t2t": {
        "1": 248387328, "2": 242696752, "3": 201105948,
        "4": 193574945, "5": 182045439, "6": 172126628,
        "7": 160567428, "8": 146259331, "9": 150617247,
        "10": 134758134, "11": 135127769, "12": 133324548,
        "13": 113566686, "14": 101161492, "15": 99753195,
        "16": 96330374, "17": 84276897, "18": 80542538,
        "19": 61707364, "20": 66210255, "21": 45090682,
        "22": 51324926, "X": 154259566, "Y": 62460029,
        "M": 16569,
    },
}

# T2T-CHM13v2.0 GenBank and RefSeq names, in chromosome order.
_T2T_GENBANK = {
    str(chrom): f"CP068{278 - chrom:03d}.2" for chrom in range(1, 23)
}
_T2T_GENBANK.update({"X": "CP068255.2", "Y": "CP086569.2"})
_T2T_REFSEQ = {
    str(chrom): f"NC_{60924 + chrom:06d}.1" for chrom in range(1, 23)
}
_T2T_REFSEQ.update({"X": "NC_060947.1", "Y": "NC_060948.1"})

T2T_ALIASES: Dict[str, str] = {}
for _aliases in (_T2T_GENBANK, _T2T_REFSEQ):
    for _chrom, _name in _aliases.items():
        T2T_ALIASES[_name] = _chrom
T2T_ALIASES["NC_012920.1"] = "M"


def normalize_assembly(value: str) -> str:
    normalized = value.lower().replace("_", "-")
    aliases = {
        "auto": "auto",
        "grch37": "grch37", "hg19": "grch37",
        "grch38": "grch38", "hg38": "grch38",
        "t2t": "t2t", "hs1": "t2t", "chm13": "t2t",
        "chm13v2.0": "t2t", "t2t-chm13v2.0": "t2t",
    }
    try:
        return aliases[normalized]
    except KeyError:
        supported = "auto, grch37/hg19, grch38/hg38, or t2t/hs1"
        raise AssemblyError(f"Unknown genome {value!r}; choose {supported}.")


def canonical_chromosome(name: str) -> str:
    if name in T2T_ALIASES:
        return T2T_ALIASES[name]
    value = name[3:] if name.lower().startswith("chr") else name
    if value.upper() in {"M", "MT"}:
        return "M"
    if value.upper() in {"X", "Y"}:
        return value.upper()
    if value.isdigit() and 1 <= int(value) <= 22:
        return str(int(value))
    return ""


def assembly_scores(references: Iterable[Tuple[str, int]]) -> Dict[str, int]:
    scores = {assembly: 0 for assembly in CHROMOSOME_LENGTHS}
    for name, length in references:
        chrom = canonical_chromosome(name)
        if not chrom:
            continue
        for assembly, lengths in CHROMOSOME_LENGTHS.items():
            if lengths.get(chrom) == length:
                scores[assembly] += 1
    return scores


def detect_assembly(references: Iterable[Tuple[str, int]]) -> str:
    scores = assembly_scores(references)
    best_score = max(scores.values())
    winners = [name for name, score in scores.items() if score == best_score]
    if best_score == 0 or len(winners) != 1:
        detail = ", ".join(f"{ASSEMBLY_LABELS[k]}={v}" for k, v in scores.items())
        raise AssemblyError(
            "Could not identify the BAM genome assembly from primary chromosome "
            f"lengths ({detail}). Pass --genome explicitly."
        )
    return winners[0]


def choose_assembly(requested: str, references: Sequence[Tuple[str, int]]) -> str:
    requested = normalize_assembly(requested)
    if requested == "auto":
        return detect_assembly(references)

    scores = assembly_scores(references)
    other_matches = [name for name, score in scores.items() if name != requested and score]
    if scores[requested] == 0 and other_matches:
        detected = max(other_matches, key=scores.get)
        raise AssemblyError(
            f"The BAM header matches {ASSEMBLY_LABELS[detected]}, not "
            f"the requested {ASSEMBLY_LABELS[requested]}."
        )
    return requested


def resolve_contig(annotation_contig: str, header_contigs: Mapping[str, int]) -> str:
    if annotation_contig in header_contigs:
        return annotation_contig
    chrom = canonical_chromosome(annotation_contig)
    for candidate in header_contigs:
        if canonical_chromosome(candidate) == chrom:
            return candidate
    raise AssemblyError(
        f"Annotation contig {annotation_contig!r} has no equivalent in the BAM header."
    )
