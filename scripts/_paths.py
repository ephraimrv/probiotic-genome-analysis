"""Repository-relative path resolution.

Every path in this repository is derived from the location of this file, never
from the working directory and never from an absolute path written into a
script. Two consequences follow, and both are the point of the arrangement:

* the repository can be moved between machines and between operating systems
  without a single edit, and
* a script produces the same result whether it is invoked from the repository
  root, from `scripts/`, or from anywhere else on the filesystem.

Paths are resolved at import time and exposed as module constants. Command-line
arguments override them; the constants are defaults, not enforcement.
"""
from __future__ import annotations

from pathlib import Path

#: Repository root — the parent of the directory holding this file.
REPO_ROOT = Path(__file__).resolve().parent.parent

DATA = REPO_ROOT / "data"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
ANNOTATION = DATA / "annotation"

RESULTS = REPO_ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
QC = RESULTS / "qc"

#: Raw data is organised per strain: data/raw/<STRAIN>/. The constants below
#: resolve to S1A, which is the pipeline's present behaviour; the --strain
#: refactor replaces them with per-strain resolution.
RAW_S1A = RAW / "S1A"

#: Assembly as delivered by SPAdes, before length filtering.
RAW_ASSEMBLY = RAW_S1A / "assembly" / "S1A_spades_contigs.fasta"

#: The 2024 Prokka annotation. Superseded by Bakta once that run exists, but
#: retained as the fallback so the pipeline is runnable from a fresh clone.
PROKKA_GENBANK = RAW_S1A / "annotation_prokka_2024-03" / "S1A_prokka.gbf"
PROKKA_CONTIG_MAP = INTERIM / "prokka_contig_map.csv"

#: Expected locations of the annotation and cluster-detection outputs.
BAKTA_GENBANK = ANNOTATION / "bakta" / "S1A.gbff"
BAKTA_PROTEINS = ANNOTATION / "bakta" / "S1A.faa"
ANTISMASH_DIR = ANNOTATION / "antismash"
BAGEL_TABLE = ANNOTATION / "bagel4" / "bagel4_results.csv"
INTERPROSCAN_TSV = ANNOTATION / "interproscan" / "S1A.faa.tsv"


def default_genbank() -> Path:
    """Bakta output when it exists, otherwise the archived Prokka annotation.

    Preferring Bakta silently is deliberate: once the re-annotation is in
    place every script should use it without a flag change. Each script
    reports the annotation it actually read, so the substitution is never
    invisible in the output.
    """
    return BAKTA_GENBANK if BAKTA_GENBANK.exists() else PROKKA_GENBANK


def rel(path: Path) -> str:
    """Render a path relative to the repository root for logging."""
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT).as_posix())
    except ValueError:
        return str(path)
