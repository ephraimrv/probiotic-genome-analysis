# `data/annotation` — annotation and cluster-detection output, not committed

The scripts in this repository locate their inputs here automatically, through
`scripts/_paths.py`. Nothing in this directory is under version control, so the
layout below has to be reproduced after a fresh clone before `run_all.sh` will
use the current annotation rather than the archived Prokka one.

| Path | Produced by | Consumed by |
|---|---|---|
| `bakta/S1A.gbff` | `scripts/run_bakta.sh` | `extract_features.py` |
| `bakta/S1A.faa` | `scripts/run_bakta.sh` | InterProScan |
| `antismash/` | antiSMASH 8, unpacked results directory | `parse_bgc.py` |
| `bagel4/bagel4_results.csv` | BAGEL4, saved by hand from the web interface | `parse_bgc.py` |
| `interproscan/S1A.faa.tsv` | InterProScan 5 | `classify_tfs.py` |

In the absence of `bakta/S1A.gbff` the pipeline falls back to the archived
Prokka annotation in `data/raw/S1A/annotation_prokka_2024-03/`. Every script prints
the annotation it read, so the fallback is recorded in the run output rather
than passing unnoticed.

The procedure for producing these files is in `docs/01_annotation_runbook.md`.
