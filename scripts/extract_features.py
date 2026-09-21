#!/usr/bin/env python3
"""Parse a bacterial GenBank annotation into tidy tables for plotting.

Annotation-agnostic: works on Prokka `.gbf`/`.gbk` and Bakta `.gbff` alike, so
swapping one annotation for the other is an input path change and nothing else.

Writes four tables:

    contigs.csv      one row per contig, with the cumulative offset used to lay
                     a draft assembly out around a circle
    features.csv     one row per annotated feature
    gc_content.csv   sliding-window GC, and its deviation from the genome mean
    gc_skew.csv      sliding-window GC skew, with cumulative skew computed
                     PER CONTIG and reset at every contig boundary

On the last point: cumulative GC skew is the standard way to locate a
replication origin, and that inference is only valid on a closed genome. This
assembly is in 57 pieces of unknown order and orientation, so a cumulative skew
running across the whole concatenation would be an artefact of contig ordering.
The reset is deliberate. Do not remove it unless the genome is closed.

Usage
-----
    python scripts/extract_features.py                  # defaults suffice
    python scripts/extract_features.py --window 2000    # or override
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

from Bio import SeqIO

# `scripts/` is not guaranteed to be on the import path: an interpreter
# started with PYTHONSAFEPATH, or -P, does not prepend the script's own
# directory. Adding it explicitly keeps the script runnable as a plain
# file from any working directory and under any interpreter setting.
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from _paths import INTERIM, PROKKA_CONTIG_MAP, TABLES, default_genbank, rel

#: Feature types carried through to the figures. Anything else is counted in
#: the summary but not written to features.csv.
KEEP_TYPES = {
    "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "misc_RNA", "regulatory",
    "repeat_region", "oriC", "oriV", "oriT", "gap", "assembly_gap",
}


def load_contig_map(path: Path | None) -> dict[str, str]:
    if path is None or not Path(path).exists():
        return {}
    with open(path) as fh:
        return {r["from_id"]: r["to_id"] for r in csv.DictReader(fh)}


def derive_contig_map(records, assembly_fasta: Path) -> dict[str, str]:
    """Match an annotation's own contig names to the stable identifiers.

    Prokka renames contigs on input, so its LOCUS names bear no relation to the
    `S1A_contig_NNN` identifiers that every table and figure in this repository
    is keyed on. The correspondence is recovered by sequence identity: each
    record is keyed on the SHA-256 of its uppercased sequence and matched
    against the renamed assembly written by `prepare_assembly.py`.

    Length was the obvious key and is the wrong one — this assembly contains
    contigs of equal length, which makes a length-keyed mapping ambiguous and
    capable of silently transposing two contigs' features. Sequence identity
    has no such failure mode, and where two contigs are byte-identical the
    choice between them cannot affect any downstream result.

    Bakta run with `--keep-contig-headers` needs none of this: its names are
    already the stable ones and the mapping comes back as the identity.
    """
    if not assembly_fasta.exists():
        return {}

    def digest(seq: str) -> str:
        return hashlib.sha256(seq.upper().encode()).hexdigest()

    stable: dict[str, str] = {}
    name, chunks = None, []
    with open(assembly_fasta) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    stable[digest("".join(chunks))] = name
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        stable[digest("".join(chunks))] = name

    mapping = {r.id: stable[digest(str(r.seq))]
               for r in records if digest(str(r.seq)) in stable}
    unmatched = [r.id for r in records if r.id not in mapping]
    if unmatched:
        raise SystemExit(
            f"cannot derive a contig mapping: {len(unmatched)} of "
            f"{len(records)} annotated contigs have no sequence match in "
            f"{assembly_fasta.name} (first: {unmatched[0]}). The annotation "
            "and the assembly are not the same object; supply an explicit "
            "--contig-map or re-run prepare_assembly.py."
        )
    return mapping


def qual(feature, key, default=""):
    v = feature.qualifiers.get(key)
    return v[0] if v else default


def gc_fraction(seq: str) -> float:
    s = seq.upper()
    g, c = s.count("G"), s.count("C")
    acgt = g + c + s.count("A") + s.count("T")
    return (g + c) / acgt if acgt else 0.0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genbank", type=Path, default=None,
                    help="GenBank input; defaults to the Bakta output when present, "
                         "otherwise the archived Prokka annotation")
    ap.add_argument("--outdir", type=Path, default=TABLES)
    ap.add_argument("--contig-map", type=Path, default=PROKKA_CONTIG_MAP,
                    help="CSV with columns from_id,to_id renaming GenBank LOCUS names")
    ap.add_argument("--window", type=int, default=5000,
                    help="sliding-window size in bp for GC tracks (default 5000)")
    ap.add_argument("--step", type=int, default=None,
                    help="window step in bp (default: equal to --window, i.e. non-overlapping)")
    ap.add_argument("--min-contig-length", type=int, default=0,
                    help="omit contigs shorter than this from the tables")
    args = ap.parse_args(argv)
    if args.genbank is None:
        args.genbank = default_genbank()

    step = args.step or args.window
    args.outdir.mkdir(parents=True, exist_ok=True)
    rename = load_contig_map(args.contig_map)

    records = [r for r in SeqIO.parse(args.genbank, "genbank")
               if len(r.seq) >= args.min_contig_length]
    if not records:
        raise SystemExit(f"no records parsed from {args.genbank}")
    if not rename:
        rename = derive_contig_map(records, INTERIM / "S1A_contigs.fasta")
    renamed = sum(1 for r in records if rename.get(r.id, r.id) != r.id)
    if not any(str(r.seq).strip("Nn") for r in records):
        raise SystemExit(f"{args.genbank} carries no sequence; GC tracks need "
                         "a GenBank file with an ORIGIN block")

    # Longest contig first, so the circle reads large-to-small clockwise.
    records.sort(key=lambda r: len(r.seq), reverse=True)

    genome_gc = gc_fraction("".join(str(r.seq) for r in records))
    total_len = sum(len(r.seq) for r in records)

    contigs, features, gc_rows, skew_rows = [], [], [], []
    offset = 0
    type_counts: dict[str, int] = {}

    for idx, rec in enumerate(records, start=1):
        cid = rename.get(rec.id, rec.id)
        seq = str(rec.seq).upper()
        n = len(seq)

        contigs.append({
            "contig_id": cid, "genbank_locus": rec.id, "order": idx,
            "length_bp": n, "cumulative_start": offset, "cumulative_end": offset + n,
            "gc_percent": round(100 * gc_fraction(seq), 3),
        })

        for feat in rec.features:
            type_counts[feat.type] = type_counts.get(feat.type, 0) + 1
            if feat.type not in KEEP_TYPES:
                continue
            start = int(feat.location.start)        # 0-based, half-open
            end = int(feat.location.end)
            strand = feat.location.strand
            features.append({
                "contig_id": cid,
                "feature_type": feat.type,
                "start": start,
                "end": end,
                "strand": "+" if strand == 1 else "-" if strand == -1 else ".",
                "length_bp": end - start,
                "locus_tag": qual(feat, "locus_tag"),
                "gene": qual(feat, "gene"),
                "product": qual(feat, "product").replace("\n", " "),
                "protein_id": qual(feat, "protein_id"),
                "ec_number": qual(feat, "EC_number"),
                "db_xref": ";".join(feat.qualifiers.get("db_xref", [])),
                "inference": ";".join(feat.qualifiers.get("inference", [])),
                "cumulative_start": offset + start,
                "cumulative_end": offset + end,
            })

        cumulative = 0.0
        for ws in range(0, max(n - args.window + 1, 1), step):
            we = min(ws + args.window, n)
            win = seq[ws:we]
            g, c = win.count("G"), win.count("C")
            gc = gc_fraction(win)
            skew = (g - c) / (g + c) if (g + c) else 0.0
            cumulative += skew
            gc_rows.append({
                "contig_id": cid, "window_start": ws, "window_end": we,
                "gc_percent": round(100 * gc, 4),
                "gc_deviation": round(100 * (gc - genome_gc), 4),
                "cumulative_start": offset + ws, "cumulative_end": offset + we,
            })
            skew_rows.append({
                "contig_id": cid, "window_start": ws, "window_end": we,
                "gc_skew": round(skew, 6),
                "cumulative_gc_skew": round(cumulative, 4),   # resets per contig
                "cumulative_start": offset + ws, "cumulative_end": offset + we,
            })
        offset += n

    def write(name, rows):
        path = args.outdir / name
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        return path

    write("contigs.csv", contigs)
    write("features.csv", features)
    write("gc_content.csv", gc_rows)
    write("gc_skew.csv", skew_rows)

    # --- sanity checks: fail loudly rather than plotting nonsense -----------
    assert sum(c["length_bp"] for c in contigs) == total_len
    assert all(f["start"] < f["end"] for f in features), "feature with start >= end"
    lens = {c["contig_id"]: c["length_bp"] for c in contigs}
    over = [f for f in features if f["end"] > lens[f["contig_id"]]]
    assert not over, f"{len(over)} features run past the end of their contig"

    kept = {}
    for f in features:
        kept[f["feature_type"]] = kept.get(f["feature_type"], 0) + 1

    print(f"genbank        : {rel(args.genbank)}")
    print(f"contig names   : {renamed}/{len(records)} mapped to stable identifiers")
    print(f"contigs        : {len(contigs)}")
    print(f"total length   : {total_len:,} bp")
    print(f"genome GC      : {100 * genome_gc:.2f} %")
    print(f"window / step  : {args.window} / {step} bp  ({len(gc_rows)} windows)")
    print(f"features kept  : {sum(kept.values())}  {kept}")
    skipped = {k: v for k, v in type_counts.items() if k not in KEEP_TYPES}
    print(f"types skipped  : {skipped}")
    print(f"wrote 4 tables to {rel(args.outdir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
