#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Assembly QC: QUAST, CheckM2 and fastANI.
#
# Run under WSL2 in the conda environment holding these tools. Three numbers
# belong in the manuscript methods whether or not a reviewer asks for them:
# assembly contiguity, genome completeness and contamination, and the ANI value
# that justifies the species name on the figure.
#
# This script must be run before run_bakta.sh, because the ANI result
# determines what species name that run may legitimately assert.
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY="${REPO_ROOT}/data/interim/S1A_contigs.fasta"
REF_DIR="${REPO_ROOT}/data/raw/reference_genomes"
OUT="${REPO_ROOT}/results/qc"
THREADS="${THREADS:-10}"

# Type-strain genomes for fastANI, downloaded from NCBI into REF_DIR as .fna.
# L. paracasei and L. casei are the pair that must be distinguished;
# L. rhamnosus serves as an outgroup sanity check.
#   GCF_000288255.1  Lacticaseibacillus paracasei subsp. paracasei JCM 8130 (type)
#   GCF_000829055.1  Lacticaseibacillus casei ATCC 393 (type)
#   GCF_000011045.1  Lacticaseibacillus rhamnosus GG

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "Filtered assembly not found at ${ASSEMBLY}."
    echo "Run:  python ${REPO_ROOT}/scripts/prepare_assembly.py"
    exit 1
fi

mkdir -p "${OUT}"

echo "==> QUAST"
quast.py "${ASSEMBLY}" -o "${OUT}/quast" --threads "${THREADS}" --labels S1A

echo "==> CheckM2"
checkm2 predict \
    --input "${ASSEMBLY}" \
    --output-directory "${OUT}/checkm2" \
    --threads "${THREADS}" \
    --force \
    -x fasta

echo "==> fastANI against type strains"
if compgen -G "${REF_DIR}/*.fna" > /dev/null; then
    ls "${REF_DIR}"/*.fna > "${OUT}/ref_list.txt"
    fastANI --query "${ASSEMBLY}" \
            --refList "${OUT}/ref_list.txt" \
            --output "${OUT}/fastani_S1A.tsv" \
            --threads "${THREADS}"
    echo "--- query, reference, ANI %, fragments mapped, total fragments ---"
    cat "${OUT}/fastani_S1A.tsv"
else
    echo "No reference genomes in ${REF_DIR}; fastANI skipped."
    echo "Download the type-strain assemblies listed above from NCBI first."
fi

echo
echo "Values to record in the manuscript:"
echo "  contigs, total length, N50, largest contig  <- results/qc/quast/report.tsv"
echo "  completeness %, contamination %             <- results/qc/checkm2/quality_report.tsv"
echo "  ANI % against each type strain              <- results/qc/fastani_S1A.tsv"
echo
echo "The species boundary is 95 % ANI. A value above it against one type strain"
echo "and clearly below against the other constitutes a defensible species call."
