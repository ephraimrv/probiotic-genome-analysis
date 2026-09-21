#!/usr/bin/env python3
"""Annotation-based prescan for bacteriocin loci, and curation worksheet builder.

This is **not** a replacement for BAGEL4 or antiSMASH. Those tools search protein
family models and precursor-peptide motifs; this script only reads the product
names an annotator already assigned. It exists for two reasons:

1.  It gives a candidate locus list immediately, from the annotation alone, that
    can be compared against the two dedicated callers when they are available.
    Agreement between three independent lines is worth more than any one.
2.  It builds the per-gene curation worksheet. Locus calls from any caller are
    boundaries; the manuscript needs gene roles, and those are assigned by hand.

Method: every CDS product is matched against a rule table of bacteriocin-locus
roles. Hits classed as `precursor`, `immunity` or `named_bacteriocin` act as
seeds; seeds within `--gap` bp of each other on the same contig merge into one
candidate locus, which is then extended by `--flank` bp. Every gene inside the
extended interval goes onto the worksheet with its predicted role and empty
columns for the curated role and the evidence behind it.

Outputs
-------
    candidate_loci.csv      one row per candidate locus
    curation_worksheet.csv  one row per gene in every candidate locus

Usage
-----
    python scripts/scan_bacteriocins.py                        # defaults suffice
    python scripts/scan_bacteriocins.py --gap 20000 --flank 8000
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

# `scripts/` is not guaranteed to be on the import path: an interpreter
# started with PYTHONSAFEPATH, or -P, does not prepend the script's own
# directory. Adding it explicitly keeps the script runnable as a plain
# file from any working directory and under any interpreter setting.
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from _paths import TABLES, rel

#: Role rules, applied in order; the first match wins. Patterns are matched
#: case-insensitively against the product string. Keep this table in the
#: repository and cite it in the methods -- it is the only thing standing
#: between a product name and a role call.
#: A bacteriocin name appears in the product string of transporters, immunity
#: proteins and regulators as often as it does in a precursor -- "Lactococcin A
#: secretion protein LcnD" is a transporter, not a lactococcin. Transport,
#: processing, immunity and regulator rules therefore run BEFORE the precursor
#: rule, and the precursor rule additionally refuses any product that carries a
#: functional word from another role. Reordering these entries will silently
#: mislabel genes.
_NOT_PRECURSOR = (r"(?!.*\b(transport|secretion|export|ATP[- ]binding|permease|"
                  r"accessory|regulator|immunity|resistance|processing|synthase|"
                  r"peptidase|protease|kinase)\b)")

ROLE_RULES: list[tuple[str, str]] = [
    ("transport", r"ATP[- ]binding.*(bacteriocin|lactococcin|lactacin|processing and transport)|"
                  r"(bacteriocin|lactococcin|lactacin|plantaricin).*(ATP[- ]binding|"
                  r"transport|secretion|export)|"
                  r"peptidase C39|C39 fam|secretion protein Lcn|HlyD|"
                  r"accessory secretion|\bComA\b|PCAT"),
    ("processing", r"peptidase.*leader|leader peptide processing|"
                   r"lanthionine|dehydratase Lan|cyclase Lan|radical SAM.*peptide"),
    ("immunity", r"immunity|CAAX amino[- ]terminal protease|\bAbi\b|abortive infection|"
                 r"self[- ]immunity"),
    ("regulator", r"response regulator|histidine kinase|accessory gene regulator|"
                  r"\bagr[ABCD]\b|quorum|autoinducer|HTH[- ]type transcriptional regulator|"
                  r"transcriptional regulator|sigma factor|\bRgg\b|LytTR"),
    ("precursor", r"bacteriocin.*(double[- ]glycine|leader)|double[- ]glycine leader|"
                  r"lanthipeptide precursor|class II[a-d]? bacteriocin|"
                  + _NOT_PRECURSOR +
                  r".*\b(plantaricin|lactacin|enterocin|pediocin|sakacin|carnobacteriocin|"
                  r"lactococcin|salivaricin|gassericin|acidocin|helveticin|nisin|"
                  r"mutacin|subtilin|bacteriocin)\b"),
    ("accessory", r"bacteriocin.*accessory|induction peptide|pheromone"),
]

#: Roles that can start a candidate locus.
SEED_ROLES = {"precursor", "immunity"}


def classify(product: str) -> str:
    p = product or ""
    for role, pattern in ROLE_RULES:
        if re.search(pattern, p, flags=re.IGNORECASE):
            return role
    if re.search(r"hypothetical|uncharacteri[sz]ed|DUF", p, flags=re.IGNORECASE):
        return "unknown"
    return "other"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", type=Path, default=TABLES / "features.csv")
    ap.add_argument("--outdir", type=Path, default=TABLES)
    ap.add_argument("--gap", type=int, default=15000,
                    help="merge seed genes closer than this into one locus (bp)")
    ap.add_argument("--flank", type=int, default=5000,
                    help="extend each locus by this much on both sides (bp)")
    ap.add_argument("--min-seeds", type=int, default=1,
                    help="discard candidate loci with fewer seed genes than this")
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.features, newline="", encoding="utf-8") as fh:
        feats = [r for r in csv.DictReader(fh) if r["feature_type"] == "CDS"]
    for f in feats:
        f["start"], f["end"] = int(f["start"]), int(f["end"])
        f["role_predicted"] = classify(f["product"])

    by_contig: dict[str, list[dict]] = {}
    for f in feats:
        by_contig.setdefault(f["contig_id"], []).append(f)
    for v in by_contig.values():
        v.sort(key=lambda f: f["start"])

    loci, worksheet = [], []
    for contig, genes in by_contig.items():
        seeds = [g for g in genes if g["role_predicted"] in SEED_ROLES]
        if not seeds:
            continue
        clusters, current = [], [seeds[0]]
        for s in seeds[1:]:
            if s["start"] - current[-1]["end"] <= args.gap:
                current.append(s)
            else:
                clusters.append(current)
                current = [s]
        clusters.append(current)

        contig_len = max(g["end"] for g in genes)
        for cl in clusters:
            if len(cl) < args.min_seeds:
                continue
            lo = max(0, min(g["start"] for g in cl) - args.flank)
            hi = min(contig_len, max(g["end"] for g in cl) + args.flank)
            inside = [g for g in genes if g["end"] > lo and g["start"] < hi]
            roles = [g["role_predicted"] for g in inside]
            locus_id = f"BGC_{len(loci) + 1:02d}"
            # A complete locus carries a precursor, an immunity gene, and either
            # a dedicated exporter or a regulator. One lone immunity-like gene
            # (a CAAX/Abi protease, say) is common across the genome and is weak
            # evidence on its own -- hence the tiering rather than a filter.
            has = {r: roles.count(r) for r in
                   ("precursor", "immunity", "transport", "processing", "regulator")}
            if has["precursor"] and has["immunity"] and (has["transport"] or has["regulator"]):
                confidence = "high"
            elif len(cl) >= 2 or (has["precursor"] and (has["transport"] or has["regulator"])):
                confidence = "medium"
            else:
                confidence = "low"
            loci.append({
                "confidence": confidence,
                "locus_id": locus_id, "contig_id": contig,
                "start": lo, "end": hi, "length_bp": hi - lo,
                "n_genes": len(inside), "n_seed_genes": len(cl),
                "n_precursor": roles.count("precursor"),
                "n_immunity": roles.count("immunity"),
                "n_transport": roles.count("transport"),
                "n_processing": roles.count("processing"),
                "n_regulator": roles.count("regulator"),
                "n_unknown": roles.count("unknown"),
                "caller": "annotation_prescan",
                "seed_products": " | ".join(sorted({g["product"] for g in cl})),
            })
            for g in inside:
                worksheet.append({
                    "locus_id": locus_id, "contig_id": contig,
                    "start": g["start"], "end": g["end"], "strand": g["strand"],
                    "locus_tag": g["locus_tag"], "gene": g["gene"],
                    "product": g["product"],
                    "length_aa": (g["end"] - g["start"]) // 3 - 1,
                    "role_predicted": g["role_predicted"],
                    "role_curated": "",      # <- fill in by hand
                    "evidence": "",          # <- BAGEL4 / antiSMASH / BLAST / Pfam hit
                    "include_in_figure": "",  # <- yes / no
                })

    for name, rows in (("candidate_loci.csv", loci), ("curation_worksheet.csv", worksheet)):
        if not rows:
            print(f"no rows for {name}")
            continue
        with open(args.outdir / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    # Plotting-ready subset. Low-confidence candidates are isolated
    # immunity-like genes, common throughout a lactobacillus genome and weak
    # evidence on their own; they remain in candidate_loci.csv but are kept off
    # the figures. Curation decisions are recorded in curation_worksheet.csv
    # and applied there, not here.
    plotted = [{"locus_id": L["locus_id"], "contig_id": L["contig_id"],
                "start": L["start"], "end": L["end"],
                "confidence": L["confidence"], "label": L["locus_id"]}
               for L in loci if L["confidence"] in ("high", "medium")]
    if plotted:
        with open(args.outdir / "bacteriocin_loci.csv", "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(plotted[0]))
            w.writeheader()
            w.writerows(plotted)

    print(f"CDS scanned      : {len(feats)}")
    counts: dict[str, int] = {}
    for f in feats:
        counts[f["role_predicted"]] = counts.get(f["role_predicted"], 0) + 1
    print(f"role calls       : { {k: v for k, v in sorted(counts.items()) if k != 'other'} }")
    print(f"candidate loci   : {len(loci)}")
    for L in loci:
        print(f"  {L['locus_id']}  {L['contig_id']}:{L['start']}-{L['end']}  "
              f"{L['n_genes']} genes  "
              f"(precursor {L['n_precursor']}, immunity {L['n_immunity']}, "
              f"transport {L['n_transport']}, regulator {L['n_regulator']})")
    print(f"worksheet rows   : {len(worksheet)}")
    print(f"wrote to {rel(args.outdir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
