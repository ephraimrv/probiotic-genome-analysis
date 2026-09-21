# Release procedure: from GitHub tag to citable DOI

Order matters in one place. Zenodo mints a DOI only for releases created after
its GitHub integration has been enabled; a release tagged before the webhook
exists does not acquire a DOI retrospectively, and recovering from that requires
cutting a second release. Step 3 therefore precedes step 4.

## 0. Before the first push

Four placeholders remain to be replaced:

| File | Field |
|---|---|
| `CITATION.cff` | `orcid`, `affiliation`, `repository-code`, `date-released` |
| `.zenodo.json` | `orcid`, `affiliation`, `related_identifiers[0].identifier` |

The staged contents should then be inspected before the first commit.

```bash
git init -b main
git add -A
git status --short
```

`.gitignore` excludes the genome data, the annotator output directories and the
reference genomes. The exclusion is deliberate: the archived record is the code
and the derived tables, not intermediates that belong in a sequence archive.
Two categories should nevertheless appear in the staged list, and are easily
lost to an over-broad ignore rule — `results/tables/*.csv`, the numbers behind
every figure, and `data/interim/uniprot_domains.json`, the cache that makes the
transcription-factor step reproducible offline.

A Zenodo DOI for analysis code does not substitute for a sequence accession.
The raw reads should be deposited in the SRA and the assembly in GenBank or ENA
before submission, and a reviewer will ask for both.

## 1. Create the repository

```bash
cd ~/python_projects/github/s1a-genome-figures
git add -A
git commit -m "Genome visualisation and regulator analysis for strain S1A"
git remote add origin git@github.com:USERNAME/s1a-genome-figures.git
git push -u origin main
```

The repository must be public before Zenodo can see it; the integration does
not list private repositories.

## 2. Check the rendered README

The three figures are referenced by relative path and render on GitHub. GitHub
also displays a "Cite this repository" control when `CITATION.cff` parses; the
absence of that control indicates a syntax error in the file.

## 3. Enable the Zenodo integration, before tagging

1. Sign in to Zenodo with the GitHub account.
2. Open the GitHub settings page in Zenodo and allow it to synchronise the repository list.
3. Locate `s1a-genome-figures` and switch the toggle on, which installs a webhook.
4. Confirm the toggle is active before proceeding.

## 4. Cut the release

```bash
git tag -a v1.0.0 -m "First archived release: figures and pipeline for strain S1A"
git push origin v1.0.0
```

A GitHub Release must then be created from that tag in the web interface. The
webhook fires on the release event rather than on the tag push, so a pushed tag
without a release produces no DOI. The release notes should state which
annotation produced the archived figures, that being the fact a later reader
most needs.

## 5. Verify the deposition, then publish

Zenodo reads `.zenodo.json` in preference to GitHub's own repository metadata.
Before publishing, the deposition page should be checked for the title and
description taken from `.zenodo.json` rather than the GitHub summary, a resolving
ORCID, the MIT licence, and an upload type of Software.

Corrections belong in `.zenodo.json`, followed by a `v1.0.1` release, rather
than in the Zenodo web form: the purpose of the file is that the metadata is
version-controlled.

## 6. Two DOIs, and which to cite

Zenodo issues both a concept DOI, which always resolves to the newest version,
and a version DOI unique to `v1.0.0`.

**The version DOI is the one that belongs in the manuscript.** A reference must
point at the exact state of the code that produced the published figures, and a
concept DOI changes meaning silently on the next release. The concept DOI is
appropriate for the README badge, where resolution to the latest version is the
intended behaviour.

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

## 7. Material for the corresponding author

Reference list entry, to be adapted to the journal's style:

> Vallente, J. E. R. (2026). *S1A genome visualisation and regulator analysis*
> (Version 1.0.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.XXXXXXX

Data availability statement:

> Analysis code, the tables underlying all genomic figures, and the figure
> source files are archived at Zenodo (https://doi.org/10.5281/zenodo.XXXXXXX).
> The draft genome assembly of strain S1A is available from [GenBank/ENA] under
> accession [ACCESSION], and the raw sequencing reads from the Sequence Read
> Archive under [ACCESSION].

Contribution statement, to be checked against the journal's CRediT categories.
The accurate description of what this repository represents is:

> [Initials]: software, formal analysis, visualization, data curation.

## 8. Revisions after submission

Bakta re-annotation, curated gene roles and reviewer requests all produce new
figures. The archived release should not be overwritten. The changes are
committed, tagged `v1.1.0`, and released again, and Zenodo mints a new version
DOI under the same concept DOI. The corresponding author must then be told which
version DOI the revised manuscript should cite; that step is easily forgotten
and yields a reference pointing at superseded figures.
