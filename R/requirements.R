# R dependencies for the bacteriocin cluster panel (Windows / RStudio).
#
# A committed renv.lock is deliberately NOT shipped from the development
# machine: a lockfile pinning conda-forge builds will not resolve against a
# CRAN-backed Windows R installation, and an unresolvable lockfile is worse
# than none. The lockfile is instead generated locally, where it then
# describes the machine on which the figures were produced:
#
#   install.packages("renv")
#   renv::init()          # scans R/ and captures what is installed
#   renv::snapshot()      # writes renv.lock -- commit that file
#
# Versions used to produce the archived figures (R 4.5.3):
#   ggplot2 4.0.3, gggenes 0.7.0, dplyr 1.2.1,
#   readr 2.2.0, ggrepel 0.9.8, scales 1.4.0

pkgs <- c("ggplot2", "gggenes", "dplyr", "readr", "ggrepel", "scales")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
invisible(lapply(pkgs, library, character.only = TRUE))
cat("R", as.character(getRversion()), "\n")
for (p in pkgs) cat(sprintf("  %-10s %s\n", p, as.character(packageVersion(p))))
