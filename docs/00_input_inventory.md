# Input inventory

Seventeen files were supplied in a single flat directory, all dated March 2024.
They constitute the complete output of one **Prokka 1.11** run executed on
Illumina BaseSpace, together with the NCBI `tbl2asn` submission files that
Prokka generates downstream. No BAGEL4, antiSMASH or Bakta output was present.

All seventeen were renamed on ingestion and sorted into `data/raw/`. The full
mapping, with SHA-256 digests and original timestamps, is in
`results/tables/rename_manifest.csv`; the resulting layout is described in
`data/raw/README.md`.

## Files consumed by the pipeline

| Purpose | File |
|---|---|
| Gene features for the circular map | `data/raw/annotation_prokka_2024-03/S1A_prokka.gbf` |
| Sequence for GC content and GC skew | the same file |
| Proteome for transcription-factor classification | `data/raw/annotation_prokka_2024-03/S1A_prokka.faa` |
| Input for Bakta re-annotation | `data/raw/assembly/S1A_spades_contigs.fasta` |

The remaining thirteen files are duplicates, NCBI submission byproducts or run
logs. They are retained for completeness and excluded from version control.

## Three observations on the duplicates

The file originally named `contigs gbf S1A.gbf` is not a second annotation. It
differs from the primary GenBank file by one line: a product name reading
`Zinc-type alcohol dehydrogenase-like protein` has been hand-edited to
`Zinc-type dehydrogenase-like protein` and wrapped mid-string. The original is
used; the edited copy is retained under `superseded/` and named to say why.

Three FASTA files describe the same assembly at different stages. The SPAdes
contigs comprise 101 sequences including fragments below 200 bp; the scaffolds
comprise 98; and the Prokka `.fna` holds the 57 contigs of at least 200 bp that
Prokka retained and renamed. The pipeline filters the unfiltered contigs itself,
so the SPAdes contig file is the input and the other two are redundant.

The `.fsa`, `.tbl`, `.sqn`, `.err` and `.val` files are `tbl2asn` submission
artefacts. They record what was prepared for NCBI and have no role in the
figures.

## Assembly characteristics

| Property | Value |
|---|---|
| Contigs retained at >= 200 bp | 57 of 101 |
| Total length | 3,061,087 bp |
| N50 / L50 | 218,690 bp / 5 |
| Largest contig | 651,111 bp |
| Ambiguous bases | 0 |
| GC content | 46.12 % |
| Median k-mer coverage | 20.3x |

Coverage analysis found no evidence of a plasmid. The thirteen contigs above
twice the median coverage are all short, 270 to 3,048 bp, with elevated GC,
consistent with collapsed rRNA operons and IS elements. Because bacteriocin loci
in lactobacilli are frequently plasmid-borne, a locus called on one of these
short high-coverage contigs would warrant particular scrutiny.

## Three properties of the supplied annotation

**No taxonomy was set.** The GenBank file records the organism as `Genus
species`, the strain as `strain` and the lineage as `Unclassified`. No
genus-specific database was therefore applied, and no product name in the file
was informed by *Lacticaseibacillus*. Prokka 1.11 also dates from 2016.
Re-annotation with Bakta is the single change that most improves everything
downstream.

**No rRNA genes were called.** The feature counts are 2,870 CDS, 53 tRNA, 1
tmRNA and 0 rRNA. Barrnap ran and found nothing. At approximately 20x coverage
SPAdes routinely collapses the repeated rRNA operons of a lactobacillus into a
single fragment or omits them, so zero is plausible for this assembly. It does
mean, however, that an earlier Proksee figure of this genome displaying an rRNA
track was either showing an empty ring or displaying an annotation other than
this one.

**Strain identity is unresolved.** That same figure was titled
*Lacticaseibacillus paracasei* BCRC-16100. BCRC 16100 is a culture-collection
accession denoting a deposited reference strain, whereas `S1A` is an in-house
isolate code. The two cannot denote the same sequenced object unless a purchased
reference strain was sequenced. The species assignment requires a stated basis
in the methods in either case: fastANI against the *L. paracasei* and *L. casei*
type-strain genomes, with the 95 % threshold as the species boundary. These two
species are difficult to separate and 16S alone does not reliably resolve them,
so ANI is the evidence a reviewer will expect.
