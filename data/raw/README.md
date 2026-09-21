# `data/raw` — source files, not committed

The genome sequence and annotation files are excluded from version control by
`.gitignore`. They are reproduced here after a fresh clone, or the scripts are
pointed elsewhere with command-line arguments.

```
assembly/
  S1A_spades_contigs.fasta      SPAdes contig assembly, 101 sequences, unfiltered
  S1A_spades_scaffolds.fasta    SPAdes scaffolds, 98 sequences; not used downstream
annotation_prokka_2024-03/
  S1A_prokka.gbf                GenBank flat file; the feature source for the figures
  S1A_prokka.fna                Prokka-renamed contigs, >= 200 bp
  S1A_prokka.faa                protein translations, 2,870 CDS
  S1A_prokka.ffn                gene nucleotide sequences
  S1A_prokka.gff                GFF3 annotation
  S1A_prokka_summary.txt        feature counts
  S1A_prokka_run.log            run log; records the Prokka version and parameters
  ncbi_submission/              tbl2asn byproducts, retained but unused
  superseded/                   a hand-edited GenBank copy with a truncated product name
```

Four files are consumed by the pipeline:

| Purpose | File |
|---|---|
| Input for `prepare_assembly.py`, and for Bakta | `assembly/S1A_spades_contigs.fasta` |
| Gene features for the circular map | `annotation_prokka_2024-03/S1A_prokka.gbf` |
| Sequence for GC content and GC skew | the same GenBank file's ORIGIN blocks |
| Proteome for transcription-factor classification | features derived from the GenBank file |

Every file was renamed on ingestion. `results/tables/rename_manifest.csv`
records the original name, the new path, the SHA-256 digest and the original
modification time, so the renaming is documented and reversible.

The archived figures under `results/figures/` were produced from the annotation
recorded in the Provenance section of the top-level `README.md`. Regenerating
them from a different annotation requires that section to be updated; a figure
whose stated input is not the one that produced it is worse than no provenance
at all.
