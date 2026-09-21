#!/usr/bin/env python3
"""Assign every predicted transcription factor in the genome to a family.

Three sources of evidence, in descending order of strength. Each protein is
tagged with the tier that actually produced its call, and the tier is carried
into the figure, so a reader can see how much of the chart rests on what.

    direct_domain   A Pfam/InterPro domain called on THIS protein -- either
                    from an annotator that records Pfam accessions, or from an
                    InterProScan run supplied with --interpro-tsv. Strongest.

    homology_domain A Pfam domain attached to the UniProt entry that the
                    annotator recorded as this protein's best database hit.
                    The domain was called on the homologue, not on the S1A
                    protein, and is transferred here by inference. Usable,
                    but not a domain call on the protein in question.

    product_name    The family is named in the annotator's product string and
                    nothing better is available. Weakest; reported, never hidden.

Proteins with regulator-like wording but no family evidence are counted as
"unassigned family" rather than dropped -- the size of that bucket is part of
the result.

Deliberately NOT counted as transcription factors, though they bind DNA:
single-stranded DNA-binding protein, nucleoid-associated proteins (HU/IHF),
transcription elongation and termination factors, DNA repair and recombination
proteins, restriction-modification systems, and topoisomerases. Sigma factors
and receiver-only response regulators are counted but kept in their own
categories, because neither is a classical sequence-specific regulator.

Usage
-----
    python scripts/classify_tfs.py             # defaults suffice
    python scripts/classify_tfs.py --offline   # cached UniProt data only
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path

# `scripts/` is not guaranteed to be on the import path: an interpreter
# started with PYTHONSAFEPATH, or -P, does not prepend the script's own
# directory. Adding it explicitly keeps the script runnable as a plain
# file from any working directory and under any interpreter setting.
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from _paths import INTERIM, INTERPROSCAN_TSV, TABLES, rel

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/accessions"
BATCH = 100

# --------------------------------------------------------------------------
# Pfam domain -> regulator family. Only DNA-binding domains and the signature
# partner domains that disambiguate them are listed. Curated for Firmicutes;
# extend it rather than loosening the product-name rules.
# --------------------------------------------------------------------------
PFAM_FAMILY: dict[str, str] = {
    "PF00126": "LysR", "PF03466": "LysR",
    "PF00440": "TetR/AcrR",
    "PF01047": "MarR", "PF12802": "MarR",
    "PF00392": "GntR", "PF07729": "GntR", "PF07701": "GntR",
    "PF00165": "AraC/XylS", "PF12833": "AraC/XylS",
    "PF00325": "Crp/Fnr", "PF13545": "Crp/Fnr",
    "PF08220": "DeoR", "PF00455": "DeoR",
    "PF00356": "LacI/GalR", "PF13377": "LacI/GalR",
    "PF03551": "PadR", "PF13601": "PadR",
    "PF01381": "XRE/Cro-CI", "PF13443": "XRE/Cro-CI", "PF12844": "XRE/Cro-CI",
    "PF04397": "LytTR",
    "PF00486": "OmpR/PhoB response regulator",
    "PF00196": "NarL/FixJ response regulator",
    "PF00072": "response regulator (receiver only)",
    "PF04545": "sigma factor", "PF04542": "sigma factor",
    "PF00140": "sigma factor", "PF04539": "sigma factor", "PF08281": "sigma factor",
    "PF01022": "ArsR/SmtB", "PF12840": "ArsR/SmtB",
    "PF13404": "Lrp/AsnC", "PF01037": "Lrp/AsnC",
    "PF00376": "MerR", "PF13411": "MerR",
    "PF01475": "Fur/Zur/PerR",
    "PF02082": "Rrf2",
    "PF01614": "IclR", "PF09339": "IclR",
    "PF01978": "TrmB",
    "PF08222": "CodY", "PF06018": "CodY",
    "PF05848": "CtsR",
    "PF01628": "HrcA",
    "PF01726": "LexA", "PF00717": "LexA",
    "PF01371": "TrpR",
    "PF03965": "BlaI/MecI",
    "PF01402": "RHH (CopG/Arc)", "PF13018": "RHH (CopG/Arc)",
    "PF08280": "Mga/PRD virulence regulator",
    "PF00874": "PRD antiterminator (BglG/MtlR)", "PF03123": "PRD antiterminator (BglG/MtlR)",
    "PF02650": "WhiA", "PF14527": "WhiA",
    "PF13412": "winged HTH, family unresolved",
    "PF13384": "winged HTH, family unresolved",
    "PF08279": "winged HTH, family unresolved",
    "PF13730": "winged HTH, family unresolved",
    "PF06971": "winged HTH, family unresolved",
    "PF04297": "winged HTH, family unresolved",
}

#: Response-regulator output domains. PF00072 on its own means a receiver
#: domain with no detected DNA-binding output -- not a transcription factor.
RR_OUTPUT = {"PF00486": "OmpR/PhoB response regulator",
             "PF00196": "NarL/FixJ response regulator",
             "PF04397": "LytTR"}

#: Product-string fallbacks, applied only when no domain evidence exists.
NAME_FAMILY: list[tuple[str, str]] = [
    (r"\bLysR\b", "LysR"), (r"\bTetR\b|\bAcrR\b", "TetR/AcrR"),
    (r"\bMarR\b", "MarR"), (r"\bGntR\b", "GntR"),
    (r"\bAraC\b|\bXylS\b", "AraC/XylS"), (r"\bCrp\b|\bFnr\b|\bCcpA\b", "Crp/Fnr"),
    (r"\bDeoR\b", "DeoR"), (r"\bLacI\b|\bGalR\b", "LacI/GalR"),
    (r"\bPadR\b", "PadR"), (r"\bXre\b|Cro/CI", "XRE/Cro-CI"),
    (r"\bLytTR\b|\bAgrA\b|accessory gene regulator protein A", "LytTR"),
    (r"\bArsR\b|\bSmtB\b", "ArsR/SmtB"), (r"\bLrp\b|\bAsnC\b", "Lrp/AsnC"),
    (r"\bMerR\b", "MerR"), (r"\bFur\b|\bZur\b|\bPerR\b", "Fur/Zur/PerR"),
    (r"\bIclR\b", "IclR"), (r"\bCodY\b", "CodY"), (r"\bCtsR\b", "CtsR"),
    (r"\bHrcA\b", "HrcA"), (r"\bLexA\b", "LexA"), (r"\bSpx\b", "Spx"),
    (r"\bNrdR\b", "NrdR"), (r"\bWhiA\b", "WhiA"), (r"\bRgg\b|\bGadR\b|\bMutR\b", "Rgg"),
    (r"\bMga\b", "Mga/PRD virulence regulator"),
    (r"\bBglG\b|\bLicT\b|\bMtlR\b|antiterminator", "PRD antiterminator (BglG/MtlR)"),
    (r"sigma[- ]?(factor|54|70)|\bECF\b", "sigma factor"),
    (r"response regulator", "response regulator (receiver only)"),
]

#: Wording that makes a protein a transcription-factor CANDIDATE.
TF_CANDIDATE = re.compile(
    r"transcription(al)?[- ]?(regulator|repressor|activator|factor)|"
    r"\bregulatory protein\b|\brepressor\b|\bactivator\b|"
    r"HTH|helix[- ]turn[- ]helix|DNA[- ]binding (protein|response)|"
    r"response regulator|sigma factor|antiterminator|\bregulon\b|"
    r"operon (regulator|repressor)|\bregulator\b", re.IGNORECASE)

#: DNA-binding but not sequence-specific regulation -- excluded outright.
NOT_TF = re.compile(
    r"single[- ]stranded DNA[- ]binding|\bssb\b|"
    r"transcription (elongation|termination|antitermination factor NusB|repair)|"
    r"\bGreA\b|\bNusA\b|\bNusB\b|\bNusG\b|\bRho\b transcription|"
    r"nucleoid[- ]associated|histone[- ]like|\bHU\b protein|integration host factor|"
    r"topoisomerase|gyrase|recombinase|resolvase|integrase|transposase|"
    r"excisionase|primase|helicase|exonuclease|endonuclease|"
    r"restriction|modification methyl|DNA polymerase|DNA ligase|"
    r"RNA polymerase (alpha|beta|omega|core)|ribosomal protein", re.IGNORECASE)


def fetch_uniprot(accessions: list[str], cache_path: Path, offline: bool) -> dict:
    cache: dict = {}
    if cache_path.exists():
        cache = json.loads(cache_path.read_text())
    todo = sorted(set(accessions) - set(cache))
    if todo and offline:
        print(f"!! --offline set; {len(todo)} accessions absent from the cache "
              "will fall back to product-name evidence")
        return cache
    if todo:
        import requests
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        for i in range(0, len(todo), BATCH):
            chunk = todo[i:i + BATCH]
            try:
                resp = requests.get(UNIPROT_URL, timeout=90, params={
                    "accessions": ",".join(chunk), "format": "json",
                    "fields": "accession,protein_name,xref_pfam,xref_interpro"})
                resp.raise_for_status()
                for entry in resp.json().get("results", []):
                    xrefs = entry.get("uniProtKBCrossReferences", [])
                    cache[entry["primaryAccession"]] = {
                        "pfam": sorted({x["id"] for x in xrefs if x["database"] == "Pfam"}),
                        "interpro": sorted({x["id"] for x in xrefs if x["database"] == "InterPro"}),
                    }
                for acc in chunk:            # remember the misses too
                    cache.setdefault(acc, {"pfam": [], "interpro": []})
            except Exception as exc:                       # noqa: BLE001
                print(f"!! UniProt batch {i // BATCH + 1} failed: "
                      f"{type(exc).__name__}: {str(exc)[:120]}", file=sys.stderr)
            print(f"   UniProt {min(i + BATCH, len(todo))}/{len(todo)}", flush=True)
            time.sleep(0.3)
        cache_path.write_text(json.dumps(cache, indent=0, sort_keys=True))
    return cache


def load_interpro(path: Path | None) -> dict[str, list[str]]:
    """InterProScan TSV -> {protein_id: [Pfam accessions]} (direct domain calls)."""
    if path is None or not path.exists():
        return {}
    out: dict[str, list[str]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 5:
                continue
            if row[3].strip().lower() in ("pfam",) and row[4].startswith("PF"):
                out.setdefault(row[0], []).append(row[4].split(".")[0])
    return {k: sorted(set(v)) for k, v in out.items()}


def family_from_pfam(pfams: list[str]) -> str:
    hits = [PFAM_FAMILY[p] for p in pfams if p in PFAM_FAMILY]
    if not hits:
        return ""
    outputs = [RR_OUTPUT[p] for p in pfams if p in RR_OUTPUT]
    if outputs:                       # an output DBD beats a bare receiver domain
        return outputs[0]
    specific = [h for h in hits if not h.startswith("winged HTH")
                and h != "response regulator (receiver only)"]
    return specific[0] if specific else hits[0]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", type=Path, default=TABLES / "features.csv")
    ap.add_argument("--outdir", type=Path, default=TABLES)
    ap.add_argument("--cache", type=Path,
                    default=INTERIM / "uniprot_domains.json")
    ap.add_argument("--interpro-tsv", type=Path, default=INTERPROSCAN_TSV)
    ap.add_argument("--loci", type=Path,
                    default=TABLES / "bacteriocin_loci.csv",
                    help="bacteriocin loci CSV; regulators inside a locus are flagged")
    ap.add_argument("--offline", action="store_true",
                    help="use only the cached UniProt data; never reach the network")
    args = ap.parse_args(argv)

    with open(args.features, newline="", encoding="utf-8") as fh:
        cds = [r for r in csv.DictReader(fh) if r["feature_type"] == "CDS"]

    direct = load_interpro(args.interpro_tsv)
    for r in cds:
        r["pfam_direct"] = sorted(set(
            re.findall(r"Pfam:(PF\d{5})", r.get("inference", ""))
            + direct.get(r.get("locus_tag", ""), [])
            + direct.get(r.get("protein_id", ""), [])))
        m = re.search(r"UniProtKB:([A-Z0-9]+)", r.get("inference", ""))
        r["uniprot"] = m.group(1) if m else ""

    cache = fetch_uniprot([r["uniprot"] for r in cds if r["uniprot"]],
                          args.cache, args.offline)

    rows = []
    for r in cds:
        product = r.get("product", "")
        if NOT_TF.search(product) or not TF_CANDIDATE.search(product):
            continue
        pf_direct = r["pfam_direct"]
        pf_hom = cache.get(r["uniprot"], {}).get("pfam", []) if r["uniprot"] else []

        family = family_from_pfam(pf_direct)
        tier, pfam_used = ("direct_domain", pf_direct) if family else ("", [])
        if not family:
            family = family_from_pfam(pf_hom)
            tier, pfam_used = ("homology_domain", pf_hom) if family else ("", [])
        if not family:
            for pattern, fam in NAME_FAMILY:
                if re.search(pattern, product, re.IGNORECASE):
                    family, tier, pfam_used = fam, "product_name", []
                    break
        if not family:
            family, tier = "unassigned family", "none"

        rows.append({
            "contig_id": r["contig_id"], "start": r["start"], "end": r["end"],
            "strand": r["strand"], "locus_tag": r["locus_tag"], "gene": r["gene"],
            "product": product, "tf_family": family, "evidence_tier": tier,
            "pfam_domains": ";".join(p for p in pfam_used if p in PFAM_FAMILY),
            "pfam_all": ";".join(pf_direct or pf_hom),
            "uniprot_hit": r["uniprot"], "in_bacteriocin_locus": "", "locus_id": "",
        })

    if args.loci and args.loci.exists():
        with open(args.loci, newline="", encoding="utf-8") as fh:
            loci = [{**L, "start": int(L["start"]), "end": int(L["end"])}
                    for L in csv.DictReader(fh)]
        for row in rows:
            for L in loci:
                if (row["contig_id"] == L["contig_id"]
                        and int(row["start"]) < L["end"] and int(row["end"]) > L["start"]):
                    row["in_bacteriocin_locus"] = "yes"
                    row["locus_id"] = L["locus_id"]
                    break

    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.outdir / "tf_assignments.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    counts: dict[tuple[str, str], int] = {}
    for r in rows:
        key = (r["tf_family"], r["evidence_tier"])
        counts[key] = counts.get(key, 0) + 1
    fam_rows = []
    for fam in sorted({k[0] for k in counts}):
        entry = {"tf_family": fam}
        for tier in ("direct_domain", "homology_domain", "product_name", "none"):
            entry[tier] = counts.get((fam, tier), 0)
        entry["total"] = sum(entry[t] for t in
                             ("direct_domain", "homology_domain", "product_name", "none"))
        entry["in_bacteriocin_locus"] = sum(
            1 for r in rows if r["tf_family"] == fam and r["in_bacteriocin_locus"] == "yes")
        fam_rows.append(entry)
    fam_rows.sort(key=lambda e: (-e["total"], e["tf_family"]))
    with open(args.outdir / "tf_family_counts.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fam_rows[0]))
        w.writeheader()
        w.writerows(fam_rows)

    tiers: dict[str, int] = {}
    for r in rows:
        tiers[r["evidence_tier"]] = tiers.get(r["evidence_tier"], 0) + 1
    unassigned = sum(1 for r in rows if r["tf_family"] == "unassigned family")
    print(f"\nCDS scanned            : {len(cds)}")
    print(f"regulator candidates   : {len(rows)}  ({100 * len(rows) / len(cds):.1f} % of CDS)")
    print(f"evidence tiers         : {tiers}")
    print(f"families assigned      : {len({r['tf_family'] for r in rows}) - (1 if unassigned else 0)}")
    print(f"unassigned family      : {unassigned}  ({100 * unassigned / len(rows):.1f} %)")
    print(f"inside bacteriocin loci: {sum(1 for r in rows if r['in_bacteriocin_locus'] == 'yes')}")
    print("\ntop families:")
    for e in fam_rows[:12]:
        print(f"  {e['total']:>4}  {e['tf_family']:<38} "
              f"direct {e['direct_domain']}, homology {e['homology_domain']}, "
              f"name {e['product_name']}")
    print(f"\nwrote {rel(args.outdir / 'tf_assignments.csv')} and "
          f"{rel(args.outdir / 'tf_family_counts.csv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
