# `data/raw` — source files, not committed

The genome sequence and annotation files are excluded from version control by
`.gitignore`. They are reproduced here after a fresh clone, or the scripts are
pointed elsewhere with command-line arguments.

Files are organised per strain. Both isolates carry the same layout, so a
script reaches either by substituting the strain directory.

```
S1A/
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
S2A/
  assembly/
    S2A_spades_contigs.fasta      SPAdes contig assembly, 1,234 sequences, unfiltered
    S2A_spades_scaffolds.fasta    SPAdes scaffolds, 1,057 sequences; not used downstream
  annotation_prokka_2024-03/
    S2A_prokka.gbf                GenBank flat file
    S2A_prokka.fna                Prokka-renamed contigs, >= 200 bp (1,057)
    S2A_prokka.faa                protein translations, 5,653 CDS
    S2A_prokka.ffn                gene nucleotide sequences
    S2A_prokka.gff                GFF3 annotation
    S2A_prokka_summary.txt        feature counts
    S2A_prokka_run.log            run log; also records a tbl2asn failure
    ncbi_submission/              tbl2asn byproducts, retained but unused
  promoters/
    S2A_promoters.fasta           nine promoter regions keyed on Prokka locus tags;
                                  generating tool and parameters unrecorded
```

Four files per strain are consumed by the pipeline, shown here for S1A:

| Purpose | File |
|---|---|
| Input for `prepare_assembly.py`, and for Bakta | `S1A/assembly/S1A_spades_contigs.fasta` |
| Gene features for the circular map | `S1A/annotation_prokka_2024-03/S1A_prokka.gbf` |
| Sequence for GC content and GC skew | the same GenBank file's ORIGIN blocks |
| Proteome for transcription-factor classification | features derived from the GenBank file |

The pipeline currently resolves these to S1A. Until the `--strain` refactor
lands, S2A is reached by passing input paths explicitly on the command line.

Every file was renamed on ingestion. `results/tables/rename_manifest.csv`
records the original name, the new path, the SHA-256 digest and the original
modification time, so the renaming is documented and reversible.

The archived figures under `results/figures/` were produced from the annotation
recorded in the Provenance section of the top-level `README.md`. Regenerating
them from a different annotation requires that section to be updated; a figure
whose stated input is not the one that produced it is worse than no provenance
at all.
