#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Linear arrow diagram of curated bacteriocin loci (gggenes).
#
# Reads the curation worksheet written by scripts/scan_bacteriocins.py. If the
# `role_curated` column has been filled in by hand it is used; otherwise the
# script falls back to `role_predicted` and says so in the subtitle, so a figure
# built on uncurated predictions can never be mistaken for a curated one.
#
# Gene fill colours are threaded to the circular map: the crimson used for
# bacteriocin loci on the whole-genome circle is the precursor colour here.
#
# Usage:
#   Rscript R/plot_clusters.R [worksheet.csv] [candidate_loci.csv] [outdir]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(gggenes); library(ggrepel)
})

# Locate the repository root by walking upwards for a marker file, so the
# script runs identically from the repository root, from R/, from an RStudio
# session opened on the .Rproj, or from a working directory on the other side
# of a filesystem share. No path below is absolute.
walk_up <- function(start, marker) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  while (dir.exists(d) && !file.exists(file.path(d, marker))) {
    parent <- dirname(d)
    if (identical(parent, d)) return(NA_character_)
    d <- parent
  }
  if (dir.exists(d)) d else NA_character_
}

# Two anchors, tried in order. Under `Rscript R/plot_clusters.R` the script's
# own path is authoritative and holds regardless of the working directory;
# inside RStudio, where --file is absent, the working directory set by the
# .Rproj serves instead.
find_repo_root <- function(marker = "CITATION.cff") {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  candidates <- character(0)
  if (length(file_arg))
    candidates <- c(candidates, dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                                      winslash = "/", mustWork = FALSE)))
  candidates <- c(candidates, getwd())
  for (start in candidates) {
    root <- walk_up(start, marker)
    if (!is.na(root)) return(root)
  }
  stop("repository root not found: no ", marker, " at or above ",
       paste(candidates, collapse = " or "))
}

ROOT      <- find_repo_root()
args      <- commandArgs(trailingOnly = TRUE)
worksheet <- if (length(args) >= 1) args[1] else file.path(ROOT, "results/tables/curation_worksheet.csv")
loci_file <- if (length(args) >= 2) args[2] else file.path(ROOT, "results/tables/candidate_loci.csv")
outdir    <- if (length(args) >= 3) args[3] else file.path(ROOT, "results/figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ROLE_COLOURS <- c(
  precursor  = "#C0392B",   # same crimson as the loci markers on the circle
  immunity   = "#E07B39",
  transport  = "#1F4E79",
  processing = "#5B9BD5",
  regulator  = "#8E6CA8",
  accessory  = "#00666E",
  unknown    = "#D9D9D9",
  other      = "#F2F2F2"
)
ROLE_LABELS <- c(
  precursor  = "precursor peptide",
  immunity   = "immunity",
  transport  = "transport / export",
  processing = "processing",
  regulator  = "regulation",
  accessory  = "accessory",
  unknown    = "hypothetical protein",
  other      = "other function"
)

ws   <- read_csv(worksheet, show_col_types = FALSE)
loci <- read_csv(loci_file, show_col_types = FALSE)

curated <- "role_curated" %in% names(ws) && any(!is.na(ws$role_curated) & ws$role_curated != "")
ws <- ws %>%
  mutate(role = if (curated) coalesce(na_if(as.character(role_curated), ""),
                                      as.character(role_predicted))
                else as.character(role_predicted))

keep <- loci %>% filter(confidence %in% c("high", "medium")) %>% pull(locus_id)
if ("include_in_figure" %in% names(ws)) {
  dropped <- ws %>% filter(tolower(as.character(include_in_figure)) == "no") %>% pull(locus_id)
  keep <- setdiff(keep, unique(dropped))
}

dat <- ws %>%
  filter(locus_id %in% keep) %>%
  group_by(locus_id) %>%
  mutate(
    origin   = min(start),
    x        = (start - origin) / 1000,
    xend     = (end   - origin) / 1000,
    forward  = strand == "+",
    role     = factor(role, levels = names(ROLE_COLOURS))
  ) %>%
  ungroup()

# Panel headers carry the locus coordinates so the kb axis is anchored.
panel_titles <- loci %>%
  filter(locus_id %in% keep) %>%
  mutate(facet = sprintf("%s  \u2014  %s : %s\u2013%s  (%s confidence)",
                         locus_id, contig_id,
                         format(start, big.mark = ","), format(end, big.mark = ","),
                         confidence)) %>%
  select(locus_id, facet)
dat <- dat %>% left_join(panel_titles, by = "locus_id")

# Label only genes with a called role or a gene symbol; 37 hypothetical
# proteins with locus tags on them would be unreadable and say nothing.
labels <- dat %>%
  filter(!role %in% c("unknown", "other") | (!is.na(gene) & gene != "")) %>%
  mutate(lab = ifelse(!is.na(gene) & gene != "", gene,
                      sub("^[A-Za-z_]*", "", locus_tag)),
         xmid = (x + xend) / 2)

p <- ggplot(dat, aes(xmin = x, xmax = xend, y = 1, fill = role, forward = forward)) +
  geom_gene_arrow(arrowhead_height = grid::unit(4.2, "mm"),
                  arrowhead_width  = grid::unit(1.6, "mm"),
                  arrow_body_height = grid::unit(3.4, "mm"),
                  colour = "grey25", linewidth = 0.22) +
  geom_text_repel(data = labels,
                  aes(x = xmid, y = 1, label = lab),
                  inherit.aes = FALSE, size = 2.3, colour = "grey15",
                  nudge_y = 0.32, direction = "x", angle = 90, hjust = 0,
                  segment.size = 0.18, segment.colour = "grey60",
                  min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  # drop = TRUE: a legend entry for a role with no genes in any panel is an
  # empty swatch a reader has to interpret. Absent roles simply do not appear.
  scale_fill_manual(values = ROLE_COLOURS, labels = ROLE_LABELS,
                    drop = TRUE, name = NULL) +
  scale_x_continuous(name = "position within locus (kb)",
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0.55, 2.1), expand = c(0, 0)) +
  facet_wrap(~facet, ncol = 1, scales = "free") +
  labs(
    # Not "predicted bacteriocin loci": the medium-confidence panels are
    # candidates under review, and some of them will not survive curation.
    title = "Candidate bacteriocin loci in the S1A draft genome",
    subtitle = if (curated)
        "Gene roles assigned by manual curation against BAGEL4 / antiSMASH evidence."
      else
        paste("Gene roles are automated predictions from annotation product",
              "names and have NOT been manually curated."),
    caption = paste("Arrows point in the direction of transcription. Genes are drawn",
                    "to scale; the kb axis restarts at each locus.")
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.title.y     = element_blank(),
    axis.text.y      = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    strip.text       = element_text(hjust = 0, face = "bold", size = 8,
                                    margin = margin(b = 3)),
    legend.position  = "bottom",
    legend.key.size  = grid::unit(3.6, "mm"),
    legend.text      = element_text(size = 7),
    plot.title       = element_text(size = 10, hjust = 0),
    plot.subtitle    = element_text(size = 7.5, colour = "grey30", hjust = 0),
    plot.caption     = element_text(size = 6.5, colour = "grey40", hjust = 0),
    plot.margin      = margin(6, 10, 6, 6)
  ) +
  guides(fill = guide_legend(nrow = 1))

h <- 1.5 + 1.35 * length(keep)
ggsave(file.path(outdir, "bacteriocin_clusters.png"), plot = p,
       width = 9.5, height = h, dpi = 600, bg = "white")
ggsave(file.path(outdir, "bacteriocin_clusters.pdf"), plot = p,
       width = 9.5, height = h, device = cairo_pdf)
# SVG alongside PDF, for consistency with the two Python figures and so that
# text remains editable in a vector editor.
ggsave(file.path(outdir, "bacteriocin_clusters.svg"), plot = p,
       width = 9.5, height = h, device = svglite::svglite)

cat(sprintf("loci plotted : %d (%s)\n", length(keep), paste(keep, collapse = ", ")))
cat(sprintf("genes drawn  : %d\n", nrow(dat)))
cat(sprintf("roles used   : %s\n", if (curated) "role_curated" else "role_predicted"))
print(table(dat$role))
cat(sprintf("wrote        : %s\n", file.path(outdir, "bacteriocin_clusters.png")))
