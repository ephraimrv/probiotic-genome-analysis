# Runbook: re-annotation and bacteriocin locus detection

Every procedure in this document runs on a local machine rather than in the
analysis session that produced this repository. Bakta requires a database not
reachable from that session, BAGEL4 exists only as a web service, and antiSMASH
installs under Linux. The plotting code was therefore built and validated
against the archived Prokka annotation, and the outputs described below drop in
as replacements without any code change.

Order matters in one place. **QC precedes Bakta**, because the fastANI result
determines what species name the Bakta command may legitimately assert.

Output locations are fixed by `scripts/_paths.py` and require no configuration:
`data/annotation/bakta/`, `data/annotation/antismash/`,
`data/annotation/bagel4/` and `data/annotation/interproscan/`. Once
`data/annotation/bakta/S1A.gbff` exists, every script uses it in preference to
the archived Prokka file and reports which annotation it read.

---

## The workflow, reviewed

The workflow originally proposed was: annotate with Prokka or Bakta; run BAGEL4
and antiSMASH to identify bacteriocin loci; curate the loci by hand, confirming
cluster boundaries and gene roles; and plot the result as a whole-genome circle
with a zoomed linear arrow diagram beside it. That is the correct workflow and
is how such papers are constructed. Four amendments apply, in order of
consequence.

**1. Bakta rather than Prokka, and for this genome the choice is not close.**
The existing annotation is Prokka 1.11 run without `--genus` or `--species`, so
no taxon-specific database was ever applied. Bakta's UniRef-derived database
yields better product names, detects small ORFs — which matters here, since the
quorum-sensing induction peptide expected in `BGC_08` is unannotated — and
writes cross-references that the transcription-factor step consumes directly.

**2. Run both locus callers, and expect disagreement.** antiSMASH matches core
biosynthetic enzymes and is strong on lanthipeptides, but can miss small
unmodified class II bacteriocins, which are what most *Lacticaseibacillus*
bacteriocins are. BAGEL4 searches precursor-peptide motifs directly. Where the
two disagree, the disagreement is the curation decision, and it is recorded
rather than resolved silently. `scripts/parse_bgc.py` merges the two with the
annotation-based prescan and flags every single-caller locus.

**3. Enable the antiSMASH TFBS finder.** It predicts transcription-factor
binding sites within cluster regions, which is *cis* evidence standing beside
the *trans* evidence of the regulator-gene counts in the transcription-factor
figure. The two together are substantially stronger than either alone.

**4. Use relaxed detection strictness.** On a 57-contig draft some clusters
necessarily sit at contig ends, and strict mode discards the partial ones.

---

## Step 1 — QC and species assignment

```bash
python scripts/prepare_assembly.py
bash scripts/run_qc.sh
```

`prepare_assembly.py` filters contigs to at least 200 bp, sorts them longest
first, assigns the stable `S1A_contig_NNN` identifiers that every table and
figure is keyed on, and records the SPAdes k-mer coverage of each in
`results/tables/contig_mapping.csv`.

`run_qc.sh` runs QUAST, CheckM2 and fastANI. The type-strain genomes for fastANI
are downloaded from NCBI into `data/raw/reference_genomes/`; the accessions are
listed in the script. The ANI value against each type strain, with 95 % as the
species boundary, is what justifies the species name used in the next step and
printed on the figure.

## Step 2 — Bakta

```bash
bash scripts/run_bakta.sh
```

The `GENUS` and `SPECIES` values in that script must be consistent with the
fastANI result. The full Bakta database is preferred over the light database:
the light database leaves substantially more proteins annotated as hypothetical,
and the transcription-factor step reads product names.

The version of Bakta and of its database should be recorded at this point; both
are printed by the script and both belong in the methods.

## Step 3 — antiSMASH

Submit `data/annotation/bakta/S1A.gbff` to the bacterial antiSMASH server, or
run a local installation. Enable KnownClusterBlast, MIBiG cluster comparison,
the TFBS finder and relaxed detection strictness. Unpack the results directory
to `data/annotation/antismash/`.

The per-region GenBank files carry the original contig coordinates, which is
what `parse_bgc.py` reads; the HTML report alone is not sufficient.

## Step 4 — BAGEL4

Submit `data/interim/S1A_contigs.fasta` to the BAGEL4 web server with all three
bacteriocin classes selected. BAGEL4 provides no API and no bulk export, so the
results table is saved by hand as
`data/annotation/bagel4/bagel4_results.csv`, with columns for contig, start,
end, class and name. The accepted column aliases are documented in
`scripts/parse_bgc.py`.

## Step 5 — reconcile and curate

```bash
python scripts/scan_bacteriocins.py
python scripts/parse_bgc.py
```

The merged table reports which callers supported each locus, the spread between
their boundary calls, and whether a locus abuts a contig edge. Loci supported by
one caller alone are retained and flagged rather than discarded, so that the
decision is explicit.

Gene roles are then curated by hand in
`results/tables/curation_worksheet.csv`. The `role_curated` and `evidence`
columns are filled in, and `include_in_figure` is set to `no` for candidates
that do not survive review. `R/plot_clusters.R` uses the curated roles once they
are present and states in the figure subtitle which of the two it used.

## Step 6 — InterProScan, optional but recommended

```bash
interproscan.sh -i data/annotation/bakta/S1A.faa -f tsv \\
    -o data/annotation/interproscan/S1A.faa.tsv -appl Pfam -goterms
python scripts/classify_tfs.py
```

This converts the majority of transcription-factor family calls from
homology-transferred domains to domains detected on the S1A proteins
themselves, which is visible in the figure as a shift in the stacked bars.

## Step 7 — regenerate

```bash
bash run_all.sh
Rscript R/plot_clusters.R
```

---

## Tools to cite

Each run prints its own version; the entries below identify the publications,
not the versions installed. DOIs should be confirmed against the journal pages
when the reference list is assembled.

| Tool | Reference |
|---|---|
| SPAdes | Bankevich et al., *J Comput Biol*, 2012 |
| QUAST | Gurevich et al., *Bioinformatics*, 2013 |
| CheckM2 | Chklovski et al., *Nature Methods*, 2023 |
| fastANI | Jain et al., *Nature Communications*, 2018 |
| Bakta | Schwengers et al., *Microbial Genomics*, 2021 |
| antiSMASH 8 | Blin et al., *Nucleic Acids Research*, 2025 |
| BAGEL4 | van Heel et al., *Nucleic Acids Research*, 2018 |
| InterProScan | Jones et al., *Bioinformatics*, 2014 |
| pyCirclize | Shimoyama, software citation |
| gggenes | Wilkins, R package |
