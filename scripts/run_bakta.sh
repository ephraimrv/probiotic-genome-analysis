#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Re-annotate the S1A draft genome with Bakta.
#
# Run under WSL2 in the conda environment holding Bakta. Output is written
# inside the repository at data/annotation/bakta/, which .gitignore excludes
# from version control, so the downstream scripts find it without configuration.
#
# The repository root is derived from this script's own location; no path below
# is absolute and nothing requires editing after the repository is moved.
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY="${REPO_ROOT}/data/interim/S1A_contigs.fasta"
OUT="${REPO_ROOT}/data/annotation/bakta"
THREADS="${THREADS:-10}"

# Database location. `bakta_db list` prints the available releases;
# `bakta_db download --output <dir> --type full` fetches one. The full database
# is preferred over the light database: the light database leaves substantially
# more proteins annotated as hypothetical, and the transcription-factor step in
# this repository reads product names.
BAKTA_DB="${BAKTA_DB:-$HOME/db/bakta}"

# --- taxonomy ---------------------------------------------------------------
# These values must be supported by the fastANI result from run_qc.sh. A
# species name that the ANI value does not justify propagates from here into
# every product name and into the figure caption.
GENUS="${GENUS:-Lacticaseibacillus}"
SPECIES="${SPECIES:-paracasei}"
STRAIN="${STRAIN:-S1A}"
LOCUS_TAG="${LOCUS_TAG:-S1A}"

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "Filtered assembly not found at ${ASSEMBLY}."
    echo "Run:  python ${REPO_ROOT}/scripts/prepare_assembly.py"
    exit 1
fi

echo "==> versions (record these for the methods section)"
bakta --version
echo "database:"
cat "${BAKTA_DB}/version.json" 2>/dev/null || echo "  no version.json in ${BAKTA_DB}"

echo
echo "==> Bakta"
mkdir -p "${OUT}"
bakta \
    --db "${BAKTA_DB}" \
    --output "${OUT}" \
    --prefix S1A \
    --genus "${GENUS}" \
    --species "${SPECIES}" \
    --strain "${STRAIN}" \
    --locus-tag "${LOCUS_TAG}" \
    --gram + \
    --min-contig-length 200 \
    --keep-contig-headers \
    --threads "${THREADS}" \
    --force \
    "${ASSEMBLY}"

echo
echo "==> feature summary"
grep -E "^(tRNA|tmRNA|rRNA|ncRNA|CDS|sORF|oriC|oriV|oriT|gap)" "${OUT}/S1A.txt" \
    || cat "${OUT}/S1A.txt"

echo
echo "Files consumed downstream, found automatically by scripts/_paths.py:"
echo "  data/annotation/bakta/S1A.gbff   circular map features"
echo "  data/annotation/bakta/S1A.fna    GC content and GC skew"
echo "  data/annotation/bakta/S1A.faa    transcription-factor classification"
echo "  data/annotation/bakta/S1A.tsv    product names with cross-references"
echo
echo "Upload to antiSMASH:  data/annotation/bakta/S1A.gbff"
echo "Upload to BAGEL4:     data/interim/S1A_contigs.fasta"
