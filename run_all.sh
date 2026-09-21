#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Regenerate every derived table and the two Python figures.
#
#   conda activate s1a-figures
#   bash run_all.sh
#
# The repository root is derived from the location of this script, so the
# command above works from any working directory and after the repository has
# been moved between machines or operating systems. No path below is absolute.
#
# Each script resolves its own defaults through scripts/_paths.py: the Bakta
# annotation under data/annotation/bakta/ when it exists, otherwise the
# archived Prokka annotation under data/raw/. Every script reports which
# annotation it read.
#
# The R cluster panel is a separate step, run on the machine that holds the R
# installation:
#   Rscript R/plot_clusters.R
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

STRAIN="${STRAIN:-S1A}"
ORGANISM="${ORGANISM:-}"          # left empty until the ANI result justifies a species

echo "repository : ${REPO_ROOT}"
echo "strain     : ${STRAIN}"
echo

echo "==> 1/6 prepare assembly"
python scripts/prepare_assembly.py --prefix "${STRAIN}"

echo "==> 2/6 extract features"
python scripts/extract_features.py

echo "==> 3/6 bacteriocin prescan"
python scripts/scan_bacteriocins.py

echo "==> 4/6 reconcile locus calls"
if [[ -d data/annotation/antismash || -f data/annotation/bagel4/bagel4_results.csv ]]; then
    python scripts/parse_bgc.py
else
    echo "    antiSMASH and BAGEL4 output absent; the prescan result stands alone."
    echo "    See docs/01_annotation_runbook.md."
fi

echo "==> 5/6 circular genome map"
python scripts/plot_circos.py --strain "${STRAIN}" ${ORGANISM:+--organism "${ORGANISM}"}

echo "==> 6/6 transcription factors"
python scripts/classify_tfs.py
python scripts/plot_tf.py --strain "${STRAIN}"

echo
echo "Complete. The cluster panel is produced separately:  Rscript R/plot_clusters.R"
