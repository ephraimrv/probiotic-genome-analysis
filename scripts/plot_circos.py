#!/usr/bin/env python3
"""Circular genome map of a draft bacterial assembly (pyCirclize).

Reads the tidy tables written by `extract_features.py` -- not a GenBank file --
so the figure is independent of which annotator produced the input.

Tracks, outermost inward:

    contig band     alternating shades, one block per contig, with a Mb axis
    CDS (+)         forward-strand coding sequences
    CDS (-)         reverse-strand coding sequences
    RNA             tRNA / rRNA / tmRNA / ncRNA
    GC content      deviation from the genome mean
    GC skew         (G-C)/(G+C), drawn per contig and reset at each boundary

Two design decisions are deliberate and should survive editing:

1.  The contig band is the outermost track and is never hidden. This assembly is
    57 contigs in arbitrary order; a reader must be able to see that adjacency on
    this circle is a drawing convention, not genomic neighbourhood.
2.  GC skew is drawn per contig with a break at every boundary. A continuous
    skew ring implies the origin/terminus inference that only a closed genome
    supports.

Usage
-----
    python scripts/plot_circos.py --strain S1A
    python scripts/plot_circos.py --strain S1A \
        --organism "Lacticaseibacillus paracasei"   # once ANI supports it
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pycirclize import Circos

# `scripts/` is not guaranteed to be on the import path: an interpreter
# started with PYTHONSAFEPATH, or -P, does not prepend the script's own
# directory. Adding it explicitly keeps the script runnable as a plain
# file from any working directory and under any interpreter setting.
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from _paths import FIGURES, TABLES, rel

# --- palette ---------------------------------------------------------------
# One hue family for coding strands, a reserved alarm hue for the bacteriocin
# loci that is used nowhere else, and a CVD-safe teal/orange pair for signed
# GC skew (never red/green).
C_CDS_FWD = "#1F4E79"
C_CDS_REV = "#5B9BD5"
C_TRNA = "#E07B39"
C_RRNA = "#6A9A3B"
C_TMRNA = "#8E6CA8"
C_NCRNA = "#B07AA1"
C_GC = "#4D4D4D"
C_SKEW_POS = "#C25E00"
C_SKEW_NEG = "#00666E"
C_LOCUS = "#C0392B"
C_CONTIG_A = "#BFBFBF"
C_CONTIG_B = "#8C8C8C"

RNA_COLOURS = {"tRNA": C_TRNA, "rRNA": C_RRNA, "tmRNA": C_TMRNA,
               "ncRNA": C_NCRNA, "misc_RNA": C_NCRNA}

# radial extents, outer -> inner
R_LOCUS = (100.5, 104.0)
R_CONTIG = (95.0, 100.0)
R_FWD = (87.5, 94.0)
R_REV = (80.0, 86.5)
R_RNA = (74.5, 79.0)
R_GC = (56.0, 72.0)
R_SKEW = (36.0, 52.0)


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def style(base: float = 8.0) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base,
        "axes.linewidth": 0.6,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def mb_formatter(v: float) -> str:
    return f"{v / 1e6:.1f} Mb"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", type=Path, default=TABLES)
    ap.add_argument("--outdir", type=Path, default=FIGURES)
    ap.add_argument("--organism", default="", help="italicised in the centre label")
    ap.add_argument("--strain", default="", help="set in roman below the organism")
    ap.add_argument("--loci", type=Path,
                    default=TABLES / "bacteriocin_loci.csv",
                    help="optional CSV: contig_id,start,end,label (bacteriocin loci)")
    ap.add_argument("--basename", default="circos_genome")
    ap.add_argument("--tick-interval", type=float, default=500_000)
    ap.add_argument("--min-label-length", type=int, default=100_000,
                    help="label contigs at least this long in the contig band")
    args = ap.parse_args(argv)

    style()
    args.outdir.mkdir(parents=True, exist_ok=True)

    contigs = read_csv(args.tables / "contigs.csv")
    features = read_csv(args.tables / "features.csv")
    gc = read_csv(args.tables / "gc_content.csv")
    skew = read_csv(args.tables / "gc_skew.csv")

    for c in contigs:
        c["length_bp"] = int(c["length_bp"])
        c["cumulative_start"] = int(c["cumulative_start"])
        c["cumulative_end"] = int(c["cumulative_end"])
    total = contigs[-1]["cumulative_end"]
    genome_gc = sum(float(c["gc_percent"]) * c["length_bp"] for c in contigs) / total

    circos = Circos(sectors={"genome": total}, start=0, end=358)
    sector = circos.sectors[0]

    # --- contig band --------------------------------------------------------
    band = sector.add_track(R_CONTIG)
    band.axis(fc="white", ec="none")
    for i, c in enumerate(contigs):
        band.rect(c["cumulative_start"], c["cumulative_end"],
                  fc=C_CONTIG_A if i % 2 == 0 else C_CONTIG_B, ec="white", lw=0.35)
    band.xticks_by_interval(args.tick_interval, outer=True, label_size=6,
                            label_formatter=mb_formatter, tick_length=1.6)
    for i, c in enumerate(contigs):
        if c["length_bp"] >= args.min_label_length:
            band.text(c["contig_id"].replace("S1A_contig_", "").lstrip("0"),
                      x=(c["cumulative_start"] + c["cumulative_end"]) / 2,
                      size=6, color="white" if i % 2 else "#1A1A1A",
                      adjust_rotation=True, orientation="vertical")

    # --- coding and RNA features -------------------------------------------
    fwd = sector.add_track(R_FWD)
    rev = sector.add_track(R_REV)
    rna = sector.add_track(R_RNA)
    for t in (fwd, rev, rna):
        t.axis(fc="#F7F7F7", ec="none")

    rna_seen: set[str] = set()
    for f in features:
        s, e = int(f["cumulative_start"]), int(f["cumulative_end"])
        if f["feature_type"] == "CDS":
            (fwd if f["strand"] == "+" else rev).rect(
                s, e, fc=C_CDS_FWD if f["strand"] == "+" else C_CDS_REV, ec="none")
        elif f["feature_type"] in RNA_COLOURS:
            rna_seen.add(f["feature_type"])
            # widen so that a 76 bp tRNA is visible at genome scale
            pad = max(0, (2500 - (e - s)) // 2)
            rna.rect(max(0, s - pad), min(total, e + pad),
                     fc=RNA_COLOURS[f["feature_type"]], ec="none")

    # --- GC content ---------------------------------------------------------
    gc_track = sector.add_track(R_GC)
    gc_track.axis(fc="none", ec="#CCCCCC", lw=0.4)
    dev = np.array([float(r["gc_deviation"]) for r in gc])
    lim = float(np.abs(dev).max())
    for c in contigs:
        rows = [r for r in gc if r["contig_id"] == c["contig_id"]]
        if len(rows) < 2:
            continue
        x = np.array([(int(r["cumulative_start"]) + int(r["cumulative_end"])) / 2 for r in rows])
        y = np.array([float(r["gc_deviation"]) for r in rows])
        gc_track.fill_between(x, y, 0, vmin=-lim, vmax=lim, fc=C_GC, ec="none", alpha=0.85)

    # --- GC skew, reset at every contig boundary ----------------------------
    sk_track = sector.add_track(R_SKEW)
    sk_track.axis(fc="none", ec="#CCCCCC", lw=0.4)
    sk_all = np.array([float(r["gc_skew"]) for r in skew])
    slim = float(np.abs(sk_all).max())
    for c in contigs:
        rows = [r for r in skew if r["contig_id"] == c["contig_id"]]
        if len(rows) < 2:
            continue
        x = np.array([(int(r["cumulative_start"]) + int(r["cumulative_end"])) / 2 for r in rows])
        y = np.array([float(r["gc_skew"]) for r in rows])
        sk_track.fill_between(x, np.clip(y, 0, None), 0, vmin=-slim, vmax=slim,
                              fc=C_SKEW_POS, ec="none")
        sk_track.fill_between(x, np.clip(y, None, 0), 0, vmin=-slim, vmax=slim,
                              fc=C_SKEW_NEG, ec="none")

    # --- bacteriocin loci ---------------------------------------------------
    loci = []
    if args.loci and args.loci.exists():
        offsets = {c["contig_id"]: c["cumulative_start"] for c in contigs}
        loci = read_csv(args.loci)
        lt = sector.add_track(R_LOCUS)
        for L in loci:
            off = offsets[L["contig_id"]]
            s, e = off + int(L["start"]), off + int(L["end"])
            lt.rect(s, max(e, s + total * 0.004), fc=C_LOCUS, ec="none")
            lt.text(L["label"], x=(s + e) / 2, r=R_LOCUS[1] + 3.0, size=6,
                    color=C_LOCUS, adjust_rotation=True, orientation="horizontal")

    # --- centre label -------------------------------------------------------
    if args.organism:
        circos.text(args.organism, r=14, size=10, style="italic", color="black")
    if args.strain:
        circos.text(f"strain {args.strain}", r=6, size=9, color="black")
    circos.text(f"{total / 1e6:.2f} Mb draft assembly\n{len(contigs)} contigs",
                r=-6, size=7, color="#555555")

    fig = circos.plotfig(dpi=100, figsize=(9, 9))

    handles = [
        Patch(fc=C_CDS_FWD, ec="none", label="CDS, forward strand"),
        Patch(fc=C_CDS_REV, ec="none", label="CDS, reverse strand"),
    ]
    handles += [Patch(fc=RNA_COLOURS[t], ec="none", label=t)
                for t in ("rRNA", "tRNA", "tmRNA", "ncRNA") if t in rna_seen]
    handles += [
        Patch(fc=C_GC, ec="none", label=f"GC content (mean {genome_gc:.1f} %)"),
        Patch(fc=C_SKEW_POS, ec="none", label="GC skew, G-rich"),
        Patch(fc=C_SKEW_NEG, ec="none", label="GC skew, C-rich"),
        Patch(fc=C_CONTIG_A, ec=C_CONTIG_B, label="contig boundaries"),
    ]
    if loci:
        handles.append(Patch(fc=C_LOCUS, ec="none", label="predicted bacteriocin locus"))
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize=7, handlelength=1.1, handleheight=1.1,
               borderaxespad=0, labelspacing=0.6)

    png = args.outdir / f"{args.basename}.png"
    fig.savefig(png, dpi=600)
    fig.savefig(args.outdir / f"{args.basename}.pdf")
    fig.savefig(args.outdir / f"{args.basename}.svg")
    print(f"contigs {len(contigs)} | total {total:,} bp | genome GC {genome_gc:.2f} %")
    print(f"GC deviation range +/-{lim:.2f} pp | GC skew range +/-{slim:.3f}")
    print(f"loci plotted: {len(loci)}")
    print(f"wrote {rel(png)} (+ .pdf, .svg)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
