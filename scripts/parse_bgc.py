#!/usr/bin/env python3
"""Reconcile bacteriocin locus calls from antiSMASH, BAGEL4 and the prescan.

Each caller answers a different question, so their outputs are merged rather
than intersected, and every merged locus records which callers supported it.
A locus found by one caller only is not discarded -- it is reported with
`n_callers = 1` so the curation decision is explicit rather than silent.

antiSMASH
    Point `--antismash` at the unpacked results directory. Region coordinates
    are read from the per-region GenBank files (`*.region*.gbk`), which carry
    the original contig coordinates in the antiSMASH structured comment.

BAGEL4
    Has no machine-readable export, so save its results table by hand as CSV
    and pass it with `--bagel`. Required columns, case-insensitive, with these
    accepted aliases:

        contig  | contig_id | scaffold | sequence
        start   | begin | from
        end     | stop | to
        class   | type | bacteriocin_class      (optional)
        name    | aoi | hit | bacteriocin       (optional)

    Coordinates are assumed 1-based inclusive, as BAGEL4 reports them, and are
    converted to the 0-based half-open convention used everywhere else here.

Prescan
    Pass `--prescan results/tables/candidate_loci.csv` to include the
    annotation-based candidates from `scan_bacteriocins.py`.

Usage
-----
    python scripts/parse_bgc.py            # reads data/annotation/{antismash,bagel4}
    python scripts/parse_bgc.py --slack 10000
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
from _paths import ANTISMASH_DIR, BAGEL_TABLE, TABLES, rel

ALIASES = {
    "contig_id": ("contig_id", "contig", "scaffold", "sequence", "seqid", "record"),
    "start": ("start", "begin", "from"),
    "end": ("end", "stop", "to"),
    "bgc_class": ("class", "type", "bacteriocin_class", "product"),
    "name": ("name", "aoi", "hit", "bacteriocin", "label"),
}


def pick(row: dict, field: str, default=""):
    lower = {k.strip().lower(): v for k, v in row.items() if k}
    for alias in ALIASES[field]:
        if alias in lower and str(lower[alias]).strip():
            return str(lower[alias]).strip()
    return default


def parse_antismash(directory: Path) -> list[dict]:
    from Bio import SeqIO

    out = []
    files = sorted(directory.rglob("*.region*.gbk"))
    for path in files:
        for rec in SeqIO.parse(path, "genbank"):
            comment = rec.annotations.get("structured_comment", {})
            data = comment.get("antiSMASH-Data", {})
            try:
                orig_start = int(str(data.get("Orig. start", "0")).replace(",", ""))
                orig_end = int(str(data.get("Orig. end", len(rec.seq))).replace(",", ""))
            except ValueError:
                orig_start, orig_end = 0, len(rec.seq)
            contig = (rec.annotations.get("original_id")
                      or data.get("Original ID", "").split()[0]
                      or rec.id)
            products, edge = [], ""
            for feat in rec.features:
                if feat.type == "region":
                    products += feat.qualifiers.get("product", [])
                    edge = (feat.qualifiers.get("contig_edge", [""])[0] or "")
            out.append({
                "caller": "antiSMASH", "contig_id": contig,
                "start": orig_start, "end": orig_end,
                "bgc_class": "/".join(sorted(set(products))),
                "name": path.stem, "contig_edge": edge,
            })
    return out


def parse_bagel(path: Path) -> list[dict]:
    delim = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh, delimiter=delim))
    out = []
    for r in rows:
        s, e = pick(r, "start"), pick(r, "end")
        if not (s and e):
            continue
        s_i, e_i = int(float(s)), int(float(e))
        if s_i > e_i:
            s_i, e_i = e_i, s_i
        out.append({
            "caller": "BAGEL4", "contig_id": pick(r, "contig_id"),
            "start": s_i - 1, "end": e_i,          # 1-based inclusive -> 0-based half-open
            "bgc_class": pick(r, "bgc_class"), "name": pick(r, "name"),
            "contig_edge": "",
        })
    return out


def parse_prescan(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return [{
            "caller": "prescan", "contig_id": r["contig_id"],
            "start": int(r["start"]), "end": int(r["end"]),
            "bgc_class": r.get("confidence", ""), "name": r["locus_id"],
            "contig_edge": "",
        } for r in csv.DictReader(fh)]


def normalise_contigs(calls: list[dict], mapping: Path | None) -> list[dict]:
    """Map caller contig names onto the stable S1A_contig_NNN identifiers."""
    if mapping is None or not mapping.exists():
        return calls
    lookup: dict[str, str] = {}
    with open(mapping, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            cid = r.get("contig_id") or r.get("to_id", "")
            for key in (cid, r.get("original_header", ""), r.get("from_id", "")):
                if key:
                    lookup[key] = cid
                    lookup[key.split()[0]] = cid
    unmapped = set()
    for c in calls:
        if c["contig_id"] in lookup:
            c["contig_id"] = lookup[c["contig_id"]]
        else:
            unmapped.add(c["contig_id"])
    if unmapped:
        print(f"!! contig names not in the mapping table, left as-is: {sorted(unmapped)}")
    return calls


def merge(calls: list[dict], slack: int) -> list[dict]:
    by_contig: dict[str, list[dict]] = {}
    for c in calls:
        by_contig.setdefault(c["contig_id"], []).append(c)

    merged = []
    for contig in sorted(by_contig):
        items = sorted(by_contig[contig], key=lambda c: (c["start"], c["end"]))
        group = [items[0]]
        for c in items[1:]:
            if c["start"] <= max(g["end"] for g in group) + slack:
                group.append(c)
            else:
                merged.append(group)
                group = [c]
        merged.append(group)

    out = []
    for i, group in enumerate(sorted(merged, key=lambda g: (g[0]["contig_id"], g[0]["start"])), 1):
        callers = sorted({g["caller"] for g in group})
        classes = sorted({g["bgc_class"] for g in group if g["bgc_class"]})
        starts = [g["start"] for g in group]
        ends = [g["end"] for g in group]
        out.append({
            "locus_id": f"LOC_{i:02d}",
            "contig_id": group[0]["contig_id"],
            "start": min(starts), "end": max(ends), "length_bp": max(ends) - min(starts),
            "callers": "+".join(callers), "n_callers": len(callers),
            "bgc_class": "; ".join(classes),
            "boundary_spread_bp": max(max(starts) - min(starts), max(ends) - min(ends)),
            "contig_edge": "yes" if any(g["contig_edge"] == "True" for g in group) else "",
            "source_calls": "; ".join(f"{g['caller']}:{g['name']}:{g['start']}-{g['end']}"
                                      for g in group),
            "label": "", "role_curation_done": "",
        })
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--antismash", type=Path, default=ANTISMASH_DIR)
    ap.add_argument("--bagel", type=Path, default=BAGEL_TABLE)
    ap.add_argument("--prescan", type=Path, default=TABLES / "candidate_loci.csv")
    ap.add_argument("--contig-mapping", type=Path,
                    default=TABLES / "contig_mapping.csv")
    ap.add_argument("--outdir", type=Path, default=TABLES)
    ap.add_argument("--slack", type=int, default=5000,
                    help="merge calls within this distance into one locus (bp)")
    args = ap.parse_args(argv)

    calls: list[dict] = []
    if args.antismash and args.antismash.exists():
        got = parse_antismash(args.antismash)
        print(f"antiSMASH regions : {len(got)}")
        calls += got
    if args.bagel and args.bagel.exists():
        got = parse_bagel(args.bagel)
        print(f"BAGEL4 areas      : {len(got)}")
        calls += got
    if args.prescan and args.prescan.exists():
        got = parse_prescan(args.prescan)
        print(f"prescan candidates: {len(got)}")
        calls += got
    if not calls:
        raise SystemExit("no caller output found -- pass at least one of "
                         "--antismash / --bagel / --prescan")

    calls = normalise_contigs(calls, args.contig_mapping)
    loci = merge(calls, args.slack)

    args.outdir.mkdir(parents=True, exist_ok=True)
    path = args.outdir / "bacteriocin_loci_merged.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(loci[0]))
        w.writeheader()
        w.writerows(loci)

    agreed = [L for L in loci if L["n_callers"] >= 2]
    print(f"\nmerged loci       : {len(loci)}  ({len(agreed)} supported by >1 caller)")
    for L in loci:
        flag = "" if L["n_callers"] > 1 else "   <- single caller, curate"
        edge = "  [contig edge]" if L["contig_edge"] else ""
        print(f"  {L['locus_id']}  {L['contig_id']}:{L['start']}-{L['end']}  "
              f"{L['callers']}  {L['bgc_class'][:40]}{edge}{flag}")
    spread = [L for L in loci if L["boundary_spread_bp"] > args.slack]
    if spread:
        print(f"\n{len(spread)} loci where callers disagree on boundaries by more "
              f"than {args.slack} bp -- resolve by hand before plotting:")
        for L in spread:
            print(f"  {L['locus_id']}: spread {L['boundary_spread_bp']} bp | {L['source_calls']}")
    print(f"\nwrote {rel(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
