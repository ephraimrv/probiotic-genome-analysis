#!/usr/bin/env python3
"""Predicted regulatory gene content by family, stacked by evidence tier.

The stacking is the point. Most family calls in a draft bacterial genome rest
on a domain transferred from a database homologue rather than one called on the
protein itself, and a bar chart that hides that difference overstates what the
genome shows. Each bar is split into the tiers defined in `classify_tfs.py`, so
the proportion of each family's count that rests on direct domain evidence is
readable off the figure.

`unassigned family` is drawn last and set apart: it is a residual, not a family,
and sorting it into the ranking would make it look like the largest family in
the genome.

Usage
-----
    python scripts/plot_tf.py --strain S1A
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# `scripts/` is not guaranteed to be on the import path: an interpreter
# started with PYTHONSAFEPATH, or -P, does not prepend the script's own
# directory. Adding it explicitly keeps the script runnable as a plain
# file from any working directory and under any interpreter setting.
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from _paths import FIGURES, TABLES, rel

TIERS = [
    ("direct_domain", "#1F4E79", "Pfam domain on this protein"),
    ("homology_domain", "#5B9BD5", "domain from best database hit"),
    ("product_name", "#BDD7EE", "family named in the annotation only"),
    ("none", "#D9D9D9", "no family evidence"),
]
C_LOCUS = "#C0392B"        # same crimson as the bacteriocin loci elsewhere
RESIDUAL = "unassigned family"


def style(base: float = 8.0) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base, "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", type=Path,
                    default=TABLES / "tf_family_counts.csv")
    ap.add_argument("--outdir", type=Path, default=FIGURES)
    ap.add_argument("--strain", default="")
    ap.add_argument("--basename", default="tf_families")
    args = ap.parse_args(argv)

    style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.counts, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("direct_domain", "homology_domain", "product_name", "none",
                  "total", "in_bacteriocin_locus"):
            r[k] = int(r[k])

    named = sorted([r for r in rows if r["tf_family"] != RESIDUAL],
                   key=lambda r: (r["total"], r["tf_family"]))
    residual = [r for r in rows if r["tf_family"] == RESIDUAL]
    ordered = residual + named               # bottom-up: residual at the bottom

    labels = [r["tf_family"] for r in ordered]
    y = list(range(len(ordered)))
    total_tf = sum(r["total"] for r in rows)
    n_assigned = total_tf - sum(r["total"] for r in residual)

    fig, ax = plt.subplots(figsize=(6.9, 0.205 * len(ordered) + 1.75))

    left = [0.0] * len(ordered)
    for key, colour, _ in TIERS:
        vals = [r[key] for r in ordered]
        ax.barh(y, vals, left=left, height=0.72, color=colour,
                edgecolor="white", linewidth=0.4, zorder=2)
        left = [a + b for a, b in zip(left, vals)]

    # Count label first, then the locus marker to its right. Placing the marker
    # on the label side of the axis puts it on top of the family name.
    widest = max(r["total"] for r in ordered)
    for yi, r in zip(y, ordered):
        ax.text(r["total"] + widest * 0.012, yi, str(r["total"]),
                va="center", ha="left", fontsize=6, color="#333333", zorder=3)
        if r["in_bacteriocin_locus"]:
            ax.plot([r["total"] + widest * 0.075], [yi], marker="o", ms=3.4,
                    color=C_LOCUS, zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.get_yticklabels()[0].set_color("#666666")     # the residual row
    ax.set_xlabel("predicted regulatory genes", fontsize=8)
    ax.set_xlim(0, max(r["total"] for r in ordered) * 1.12)
    ax.set_ylim(-0.8, len(ordered) - 0.3)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, color="#EDEDED", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    strain = f" in {args.strain}" if args.strain else ""
    ax.set_title(
        f"{n_assigned} of {total_tf} predicted regulators{strain} fall into "
        f"{len(named)} families",
        fontsize=9, loc="left", pad=16)
    fig.text(0.0, 1.0,
             "Predicted regulatory gene content from genome annotation. "
             "This is gene content, not transcription-factor activity:\n"
             "no expression data supports these counts, and family calls "
             "rest mostly on homology rather than direct domain detection.",
             fontsize=6.5, color="#555555", va="bottom", ha="left",
             transform=ax.transAxes)

    handles = [mpl.patches.Patch(fc=c, ec="white", lw=0.4, label=lab)
               for _, c, lab in TIERS]
    handles.append(mpl.lines.Line2D([], [], marker="o", ls="none", ms=3.6,
                                    color=C_LOCUS,
                                    label="family present in a candidate bacteriocin locus"))
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=6.2,
              handlelength=1.1, borderaxespad=0.4, labelspacing=0.45)

    fig.savefig(args.outdir / f"{args.basename}.png", dpi=600)
    fig.savefig(args.outdir / f"{args.basename}.pdf")
    fig.savefig(args.outdir / f"{args.basename}.svg")

    print(f"families plotted : {len(named)} (+ residual)")
    print(f"total regulators : {total_tf}  | assigned to a family: {n_assigned}")
    for key, _, lab in TIERS:
        print(f"  {sum(r[key] for r in rows):>4}  {lab}")
    print(f"families with a locus regulator: "
          f"{[r['tf_family'] for r in rows if r['in_bacteriocin_locus']]}")
    print(f"wrote {rel(args.outdir / (args.basename + '.png'))} (+ .pdf, .svg)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
