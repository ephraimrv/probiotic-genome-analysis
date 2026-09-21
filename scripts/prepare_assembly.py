#!/usr/bin/env python3
"""Prepare the raw SPAdes assembly for annotation.

Filters contigs by length, sorts them longest-first, gives them stable short
names, and records a mapping table back to the original SPAdes headers (which
carry the k-mer coverage values the QC section of the manuscript needs).

Stable contig names matter downstream: every figure, table and locus coordinate
in this repository is keyed on them, so they must not change between runs.

Standard library only -- runs anywhere, no environment needed.

Usage
-----
    python scripts/prepare_assembly.py            # defaults suffice
    python scripts/prepare_assembly.py --min-length 500   # or override
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# `scripts/` is not guaranteed to be on the import path: an interpreter
# started with PYTHONSAFEPATH, or -P, does not prepend the script's own
# directory. Adding it explicitly keeps the script runnable as a plain
# file from any working directory and under any interpreter setting.
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from _paths import INTERIM, RAW_ASSEMBLY, TABLES, rel

SPADES_HEADER = re.compile(
    r"^NODE_(?P<node>\d+)_length_(?P<length>\d+)_cov_(?P<cov>[\d.]+)"
)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def wrap(seq: str, width: int = 60):
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, default=RAW_ASSEMBLY)
    ap.add_argument("--output", type=Path, default=INTERIM / "S1A_contigs.fasta")
    ap.add_argument("--mapping", type=Path, default=TABLES / "contig_mapping.csv")
    ap.add_argument("--prefix", default="S1A")
    ap.add_argument("--min-length", type=int, default=200,
                    help="drop contigs shorter than this (default 200, the NCBI floor)")
    args = ap.parse_args(argv)

    records = read_fasta(args.input)
    if not records:
        sys.exit(f"no sequences read from {args.input}")

    kept = [(h, s) for h, s in records if len(s) >= args.min_length]
    kept.sort(key=lambda hs: len(hs[1]), reverse=True)
    dropped = len(records) - len(kept)
    width = max(3, len(str(len(kept))))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.mapping.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(args.output, "w", newline="\n") as out:
        for i, (header, seq) in enumerate(kept, start=1):
            new_id = f"{args.prefix}_contig_{i:0{width}d}"
            m = SPADES_HEADER.match(header)
            gc = sum(seq.upper().count(b) for b in "GC")
            rows.append({
                "contig_id": new_id,
                "original_header": header,
                "spades_node": m.group("node") if m else "",
                "length_bp": len(seq),
                "kmer_coverage": m.group("cov") if m else "",
                "gc_percent": round(100 * gc / len(seq), 2),
                "n_count": seq.upper().count("N"),
            })
            out.write(f">{new_id}\n")
            for chunk in wrap(seq):
                out.write(chunk + "\n")

    with open(args.mapping, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    total = sum(r["length_bp"] for r in rows)
    lengths = sorted((r["length_bp"] for r in rows), reverse=True)
    cum, n50, l50 = 0, lengths[-1], len(lengths)
    for i, x in enumerate(lengths, start=1):
        cum += x
        if cum >= total / 2:
            n50, l50 = x, i
            break
    gc_total = sum(r["gc_percent"] * r["length_bp"] for r in rows) / total

    print(f"input          : {rel(args.input)}")
    print(f"sequences read : {len(records)}")
    print(f"dropped (<{args.min_length} bp): {dropped}")
    print(f"contigs kept   : {len(rows)}")
    print(f"total length   : {total:,} bp")
    print(f"N50 / L50      : {n50:,} bp / {l50}")
    print(f"largest contig : {lengths[0]:,} bp")
    print(f"GC             : {gc_total:.2f} %")
    print(f"ambiguous bases: {sum(r['n_count'] for r in rows)}")
    print(f"wrote          : {rel(args.output)}")
    print(f"wrote          : {rel(args.mapping)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
